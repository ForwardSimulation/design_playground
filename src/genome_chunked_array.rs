use crate::common::Mutation;

// Mutations will be stored in blocks of 64
// indexes
const CHUNK_SIZE: usize = 64;

struct Chunk {
    start: u32,
    stop: u32,
}

struct MutationChunks {
    mutations: Vec<u32>, // indexes into mutation vector.
    // u32::MAX is treated as a NULL "sentinel"
    chunks: Vec<Chunk>, // Is this needed?
}

impl MutationChunks {
    const CHUNK_SIZE: usize = CHUNK_SIZE; // Maybe this should be here?
}

struct Genomes {
    chunks: Vec<u32>,   // each genome is a range of chunk indexes
    starts: Vec<usize>, // For each genome, where is its first chunk?
    stops: Vec<usize>,  // One past last chunk such that a genome is chunks[starts[i]..stops[i]]
}

// What we really need:
// 1. Do we need a new chunk?
// 2. If yes, push (or recycle) one.
// 3. If no, fill what we can of the current (last?) chunk.
fn add_mutations(mutation_keys: &[u32], genome: usize, genomes: &mut Genomes) {
    genomes.starts.push(0);
    genomes.stops.push(CHUNK_SIZE);
    genomes.chunks.extend_from_slice(mutation_keys);
    for _ in 0..CHUNK_SIZE - mutation_keys.len() {
        genomes.chunks.push(u32::MAX);
    }
}

impl Genomes {}

#[cfg(test)]
mod development_tests {
    use super::*;

    #[test]
    fn test_add_mutations() {
        let mut genomes = Genomes {
            chunks: vec![],
            starts: vec![],
            stops: vec![],
        };

        let mutations = [1, 2, 3];
        add_mutations(&mutations, 0, &mut genomes);
        assert_eq!(genomes.starts.len(), 1);
        assert_eq!(genomes.stops.len(), 1);
        assert_eq!(genomes.stops[0] - genomes.starts[0], CHUNK_SIZE);
        assert_eq!(genomes.chunks.len(), CHUNK_SIZE);
    }
}

#[cfg(test)]
mod sinful_tests {
    use std::num::NonZeroU32;

    use super::*;

    #[test]
    fn test_sizes() {
        assert_eq!(
            std::mem::size_of::<Chunk>(),
            std::mem::size_of::<rclite::Rc<Chunk>>()
        );
        assert_eq!(
            std::mem::size_of::<Chunk>(),
            std::mem::size_of::<Option<rclite::Rc<Chunk>>>()
        );
        assert_eq!(
            std::mem::size_of::<Option<NonZeroU32>>(),
            std::mem::size_of::<NonZeroU32>()
        );
    }

    #[test]
    fn test_non_zero_int_types() {
        let x = NonZeroU32::new(1).unwrap();
        let y = x.checked_add(1).unwrap();
        assert_eq!(y, 2.try_into().unwrap());
    }
}
