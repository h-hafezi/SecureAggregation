use ark_ec::AffineRepr;
use ark_ff::{PrimeField, Zero};
use ark_std::vec::Vec;
use ark_ec::pairing::Pairing;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator};

use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use crate::merkle_tree::mimc::MiMC;
use crate::ring::cast_field;

/// Array-based Merkle tree using Poseidon hash outputs in field `F`.
#[derive(Clone, Debug)]
pub struct MerkleTree<E>
where
    E: Pairing,
{
    /// number of leaves (must be power of two)
    pub leaf_count: usize,
    /// tree array (1-indexed). length == 2 * leaf_count
    /// valid indices: 1 .. 2*leaf_count-1
    pub tree: Vec<E::ScalarField>,
}

const ROUNDS: usize = 91;

#[derive(Clone, Debug)]
pub struct MerkleProof<E>
where
    E: Pairing,
{
    /// sibling path from leaf level up to (but not including) root.
    /// siblings[0] is the sibling at leaf level; siblings.last() is the sibling at the top level.
    pub siblings: Vec<E::ScalarField>,
    /// index of the leaf (0-based)
    pub index: usize,
}

impl<E> MerkleTree<E>
where
    E: Pairing,
{
    /// Build a Merkle tree from `leaves` (Vec<G>).
    /// `leaves.len()` must be a power of two.
    pub fn build(leaves: Vec<E::G1Affine>) -> Self
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: ark_ff::PrimeField,
        <E as Pairing>::ScalarField: ark_crypto_primitives::sponge::Absorb
    {
        assert!(!leaves.is_empty(), "leaves must not be empty");
        let leaf_count = leaves.len();
        assert!(leaf_count.is_power_of_two(), "leaf count must be power of two");

        // Tree vector uses 1-indexing; index 0 unused.
        let mut tree: Vec<E::ScalarField> = vec![E::ScalarField::zero(); 2 * leaf_count];

        // Hash leaves into field elements and put at positions leaf_count .. 2*leaf_count-1
        for (i, g) in leaves.iter().enumerate() {
            let leaf_hash = Self::hash_leaf(g);
            tree[leaf_count + i] = leaf_hash;
        }

        // Build internal nodes
        for idx in (1..leaf_count).rev() {
            let left = tree[2 * idx];
            let right = tree[2 * idx + 1];
            let node_hash = Self::hash_pair(&left, &right);
            tree[idx] = node_hash;
        }

        MerkleTree {
            leaf_count,
            tree,
        }
    }

    /// Return the root (field element)
    pub fn root(&self) -> E::ScalarField {
        self.tree[1]
    }

    /// Generate inclusion proof for leaf at `index` (0-based).
    /// Returns MerkleProof with the sibling path (bottom-up).
    pub fn gen_proof(&self, index: usize) -> MerkleProof<E> {
        assert!(index < self.leaf_count, "leaf index out of bounds");
        let mut siblings = Vec::with_capacity(self.leaf_count.trailing_zeros() as usize);
        let mut pos = self.leaf_count + index; // 1-indexed position in tree vector

        while pos > 1 {
            let sibling_pos = if pos % 2 == 0 { pos + 1 } else { pos - 1 };
            siblings.push(self.tree[sibling_pos]);
            pos /= 2; // move to parent
        }

        MerkleProof { siblings, index }
    }

    /// Verify a proof for a given leaf and expected root.
    /// Returns true if valid.
    pub fn verify_proof(leaf: &E::G1Affine, proof: &MerkleProof<E>, expected_root: E::ScalarField) -> bool
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: ark_ff::PrimeField,
        <E as Pairing>::ScalarField: ark_crypto_primitives::sponge::Absorb
    {
        // Recompute leaf hash
        let mut cur = Self::hash_leaf(leaf);
        let mut idx = proof.index;

        for sibling in proof.siblings.iter() {
            if idx % 2 == 0 {
                // leaf is left child
                cur = Self::hash_pair(&cur, sibling);
            } else {
                // leaf is right child
                cur = Self::hash_pair(sibling, &cur);
            }
            idx /= 2;
        }

        cur == expected_root
    }

    /// Hash a G element into a Poseidon-field hash F (leaf hashing)
    fn hash_leaf(g: &E::G1Affine) -> E::ScalarField
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
        <E as Pairing>::ScalarField: ark_crypto_primitives::sponge::Absorb
    {
        let (x, y) = g.xy().unwrap();
        let x_native = cast_field::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(x);
        let y_native = cast_field::<<<E as Pairing>::G1Affine as AffineRepr>::BaseField, E::ScalarField>(y);

        // Apply Poseidon sponge to get the final leaf hash. We create a fresh PoseidonHash.
        let mimc = MiMC::<E::ScalarField>::new(ROUNDS);
        mimc.hash_elements(&[x_native, y_native])
    }

    /// Hash two field elements (left, right) into parent node hash with Poseidon
    fn hash_pair(left: &E::ScalarField, right: &E::ScalarField) -> E::ScalarField
    where
        <E as Pairing>::ScalarField: ark_crypto_primitives::sponge::Absorb
    {
        let mimc = MiMC::<E::ScalarField>::new(ROUNDS);
        mimc.hash_elements(&[*left, *right])
    }
}

impl<E> MerkleTree<E>
where
    E: Pairing + Sync + Send,
{
    /// Parallelized build: split `leaves` into `num_subtrees` equal chunks (num_subtrees must be power of two
    /// and divide the leaf_count), compute each chunk's leaf-hashes in parallel, then assemble the full tree.
    pub fn build_parallel(
        leaves: Vec<E::G1Affine>,
        num_subtrees: usize,
    ) -> Self
    where
        <<E as Pairing>::G1Affine as AffineRepr>::BaseField: PrimeField,
        <E as Pairing>::ScalarField: ark_crypto_primitives::sponge::Absorb + Send + Sync + Clone,
        E::G1Affine: Send + Sync + Clone,
    {
        assert!(!leaves.is_empty(), "leaves must not be empty");
        let leaf_count = leaves.len();
        assert!(leaf_count.is_power_of_two(), "leaf count must be power of two");
        assert!(num_subtrees.is_power_of_two(), "num_subtrees must be power of two");
        assert_eq!(leaf_count % num_subtrees, 0, "num_subtrees must evenly divide leaf_count");

        let chunk_size = leaf_count / num_subtrees;

        // Partition leaves into chunks of size chunk_size
        let chunks: Vec<Vec<E::G1Affine>> = leaves
            .chunks(chunk_size)
            .map(|c| c.to_vec())
            .collect();

        // Build each subtree in parallel using the existing `build` function.
        // Each subtree is a full MerkleTree for its chunk.
        let subtrees: Vec<MerkleTree<E>> = chunks
            .into_par_iter()
            .map(|chunk_leaves| MerkleTree::build(chunk_leaves))
            .collect();

        // Flatten leaf hashes from subtrees into the final leaves vector (ordered)
        let mut flat_leaf_hashes: Vec<E::ScalarField> = Vec::with_capacity(leaf_count);
        for subtree in subtrees.iter() {
            // subtree.leaf_count == chunk_size
            let sc = subtree.leaf_count;
            // leaves are stored at indices sc .. 2*sc - 1 in the subtree.tree
            flat_leaf_hashes.extend_from_slice(&subtree.tree[sc..2 * sc]);
        }
        assert_eq!(flat_leaf_hashes.len(), leaf_count);

        // Build final tree array: 1-indexed, size 2*leaf_count
        let mut tree: Vec<E::ScalarField> = vec![E::ScalarField::zero(); 2 * leaf_count];

        // Put leaf hashes at positions leaf_count .. 2*leaf_count-1
        for (i, leaf_hash) in flat_leaf_hashes.into_iter().enumerate() {
            tree[leaf_count + i] = leaf_hash;
        }

        // Compute internal nodes bottom-up (sequential). This uses the leaves we filled above.
        for idx in (1..leaf_count).rev() {
            let left = tree[2 * idx];
            let right = tree[2 * idx + 1];
            let node_hash = Self::hash_pair(&left, &right);
            tree[idx] = node_hash;
        }

        MerkleTree {
            leaf_count,
            tree,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{Bls12_381, G1Affine, G1Projective, Fr};
    use ark_ec::CurveGroup;
    use ark_std::UniformRand;
    use rand::thread_rng;

    #[test]
    fn test_merkle_tree_poseidon() {
        let mut rng = thread_rng();
        const N: usize = 32; // must be power of two

        // Generate random G1Affine leaves
        let leaves: Vec<G1Affine> = (0..N)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect();

        // Build the tree
        let tree = MerkleTree::<Bls12_381>::build_parallel(leaves.clone(), 4);
        let root = tree.root();
        assert!(!root.is_zero(), "root should not be zero");

        // Verify inclusion for each leaf
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.gen_proof(i);
            assert!(
                MerkleTree::<Bls12_381>::verify_proof(leaf, &proof, root),
                "proof should verify for leaf {}",
                i
            );
        }

        // Negative test: tamper one sibling in a valid proof
        let mut bad_proof = tree.gen_proof(0);
        bad_proof.siblings[0] += Fr::from(1u64); // corrupt sibling
        assert!(
            !MerkleTree::<Bls12_381>::verify_proof(&leaves[0], &bad_proof, root),
            "tampered proof must fail"
        );
    }
}
