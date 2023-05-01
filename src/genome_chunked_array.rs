use crate::common::Mutation;

// Mutations will be stored in blocks of 64
// indexes
const CHUNK_SIZE: usize = 64;

struct Chunk {
    start: u32,
    stop: u32,
}

#[derive(Default)]
struct MutationChunks {
    mutations: Vec<u32>, // indexes into mutation vector.
    // u32::MAX is treated as a NULL "sentinel"
    chunks: Vec<Chunk>, // Is this needed?
    counts: Vec<u32>,
    queue: Vec<usize>,
}

impl MutationChunks {
    const CHUNK_SIZE: usize = CHUNK_SIZE; // Maybe this should be here?

    fn new_chunk_mut(&mut self) -> (usize, &mut [u32]) {
        assert_eq!(self.mutations.len() / Self::CHUNK_SIZE, 0);

        match self.queue.pop() {
            Some(index) => {
                let mslice =
                    &mut self.mutations[index * CHUNK_SIZE..index * CHUNK_SIZE + CHUNK_SIZE];
                mslice.fill(u32::MAX);
                (index, mslice)
            }
            None => {
                let id = self.mutations.len() / Self::CHUNK_SIZE;
                self.mutations
                    .resize(self.mutations.len() + Self::CHUNK_SIZE, u32::MAX);
                (
                    id,
                    &mut self.mutations[id * CHUNK_SIZE..id * CHUNK_SIZE + CHUNK_SIZE],
                )
            }
        }
    }
}

struct Genomes {
    chunks: Vec<u32>,   // each genome is a range of chunk indexes
    starts: Vec<usize>, // For each genome, where is its first chunk?
    stops: Vec<usize>,  // One past last chunk such that a genome is chunks[starts[i]..stops[i]]
}

impl Genomes {
    fn is_empty(&self) -> bool {
        assert_eq!(self.starts.len(), self.stops.len());
        self.starts.is_empty()
    }
}

// What we really need:
// 1. Do we need a new chunk?
// 2. If yes, push (or recycle) one.
// 3. If no, fill what we can of the current (last?) chunk.
// Actually, a lot of this is wrong:
// 1. We don't know what a "genome" is, yet, in memory.
// 2. So, are we adding to a genome, updating an existing
//    genome, or what?
// So what we actually really need is:
// 1. To know the parental genome identifiers
// 2. To be able to get new chunks, which means that MutationChunks
//    has to be part of the conversation.
// 3. Need partition searching to copy parental stuff over.
// 4. Etc..
fn generate_offspring_genome(
    mutation_keys: &[u32], // The indexes of the new mutations
    genome: usize,         // The index of the genome that will "get" the new mutations
    genomes: &mut Genomes, // The output
) {
    if genomes.is_empty() {
        genomes.starts.push(0);
        genomes.stops.push(CHUNK_SIZE);
        genomes.chunks.extend_from_slice(mutation_keys);
        for _ in 0..CHUNK_SIZE - mutation_keys.len() {
            genomes.chunks.push(u32::MAX);
        }
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
        generate_offspring_genome(&mutations, 0, &mut genomes);
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
