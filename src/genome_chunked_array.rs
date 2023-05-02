use crate::common::Mutation;

// Mutations will be stored in blocks of 64
// indexes
const CHUNK_SIZE: usize = 64;

// Pseudocode of what an OO
// version may look like
mod pseudocode {
    struct Genome {
        chunks: Vec<u32>, // indexes into DiploidPopulation::mutation_chunks::chunks
        count: u32,
    }

    struct MutationChunk {
        mutations: [u32; super::CHUNK_SIZE], // indexes into DiploidPopulation::mutations
        count: u32,
    }

    impl MutationChunk {
        fn new_empty() -> Self {
            Self {
                mutations: [u32::MAX; super::CHUNK_SIZE],
                count: 0,
            }
        }
    }

    struct MutationChunks {
        chunks: Vec<MutationChunk>,
    }

    struct DiploidPopulation {
        genomes: Vec<Genome>,
        mutation_chunks: MutationChunks,
        mutations: Vec<super::Mutation>,
    }
}

struct Chunk {
    start: u32,
    stop: u32,
}

#[derive(Default)]
pub struct MutationChunks {
    mutations: Vec<u32>, // indexes into mutation vector.
    // u32::MAX is treated as a NULL "sentinel"
    chunks: Vec<Chunk>, // Is this needed?
    counts: Vec<u32>,
    queue: Vec<usize>,
}

impl MutationChunks {
    const CHUNK_SIZE: usize = CHUNK_SIZE; // Maybe this should be here?

    pub fn chunk(&self, at: usize) -> &[u32; CHUNK_SIZE] {
        let s = &self.mutations[at * CHUNK_SIZE..at * CHUNK_SIZE + CHUNK_SIZE];
        s.try_into().unwrap()
    }

    pub fn chunk_mut(&mut self, at: usize) -> &mut [u32; CHUNK_SIZE] {
        let s = &mut self.mutations[at * CHUNK_SIZE..at * CHUNK_SIZE + CHUNK_SIZE];
        s.try_into().unwrap()
    }

    // The None path may not be the most efficient
    fn new_chunk_mut(&mut self) -> (usize, &mut [u32; CHUNK_SIZE]) {
        assert_eq!(self.mutations.len() / Self::CHUNK_SIZE, 0);

        match self.queue.pop() {
            Some(index) => (index, self.chunk_mut(index)),
            None => {
                let id = self.mutations.len() / Self::CHUNK_SIZE;
                self.mutations
                    .resize(self.mutations.len() + Self::CHUNK_SIZE, u32::MAX);
                (id, self.chunk_mut(id))
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

    fn genome(&self, index: usize) -> &[u32] {
        &self.chunks[index * CHUNK_SIZE..index * CHUNK_SIZE + CHUNK_SIZE]
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
    // NOTE: we are already bowing
    // to the borrow checker:
    // we cannot have parental/offspring
    // genomes in the same container.
    // Do "chunks" and "genomes" need to
    // be in the same struct?
    parents: (&[u32], &[u32]),            // parental genomes
    mutation_chunks: &mut MutationChunks, // Output chunks
    genomes: &mut Genomes,                // Output genomes
) {
    let parent_one_genome = parents.0;
    let parent_two_genome = parents.1;
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

        let mut chunks = MutationChunks::default();
        let mutations = [1, 2, 3];
        generate_offspring_genome(&mutations, (&[], &[]), &mut chunks, &mut genomes);
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

    //#[test]
    //fn test_arrayref() {
    //    fn mod_array(x: &mut [i32; 3]) {
    //        x[2] = 4;
    //    }
    //    use arrayref::array_mut_ref;
    //    let mut x = vec![1, 2, 3];
    //    mod_array(array_mut_ref![x, 0, 3]);
    //    assert_eq!(x[2], 4);
    //}

    #[test]
    fn test_slice_try_into_array() {
        let mut v = vec![1; 64];
        let vs = &mut v[0..32];
        let s: &mut [i32; 32] = vs.try_into().unwrap();
        s[10] = 2;
        assert_eq!(v[10], 2);
        let vs = &mut v[0..64];
        let s: &mut [i32; 64] = vs.try_into().unwrap();
        s[10] = 3;
        assert_eq!(v[10], 3);

        let mut v = vec![1; 256];
        let vs = &mut v[9..9 + 64];
        let s: &mut [i32; 64] = vs.try_into().unwrap();
        s[1] = 3;
        assert_eq!(v[10], 3);
    }

    #[test]
    fn test_copy_free_array_from_struct() {
        struct X {
            x: Vec<i32>,
        }

        impl X {
            fn chunk(&mut self) -> &mut [i32; 4] {
                let s = &mut self.x[0..4];
                s.try_into().unwrap()
            }
        }

        let mut x = X {
            x: vec![0, 1, 2, 3],
        };
        let c = x.chunk();
        c[0] += 10;
        assert_eq!(x.x[0], 10);
    }
}
