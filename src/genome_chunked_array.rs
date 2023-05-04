//! This is a first draft of storing genomes
//! in "chunks" of a fixed length. A "chunk"
//! is a number of mutations.
//!
//! For this first pass, we will only support
//! efficient operations for non-overlapping
//! generations. Like genome_array.rs, each
//! new generation of genomes will be a new
//! container. Also like genome_array, this
//! design (probably) has efficiency issues
//! for overlapping generations. We will deal
//! with that later.

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

#[derive(Default)]
struct Mutations {
    mutations: Vec<Mutation>,
    counts: Vec<u32>,
}

#[derive(Default)]
struct MutationChunks {
    mutation_ids: Vec<u32>, // indexes into mutation vector.
    counts: Vec<u32>,
    // How many valid mutation ids are in a chunk.
    // The rest must be the sentinel value u32::MAX
    // NOTE: defining this as "how many" is akin to
    // the len() of a Vec.  That definition may
    // make doing things like getting
    // the last position out of a chunk awkward.
    occupancy: Vec<i32>,
    // indexes where chunks have count of 0.
    // We can use these to "reallocate" new chunks.
    queue: Vec<usize>,
}

impl MutationChunks {
    fn new_chunk(&mut self) -> usize {
        match self.queue.pop() {
            Some(index) => {
                assert_eq!(self.counts[index], 0);
                self.mutation_ids[index * CHUNK_SIZE..index * CHUNK_SIZE + CHUNK_SIZE]
                    .fill(u32::MAX);
                self.occupancy[index] = 0;
                index
            }
            None => {
                let index = self.mutation_ids.len();
                self.mutation_ids.resize(index + CHUNK_SIZE, u32::MAX);
                self.occupancy.push(0);
                self.counts.push(0);
                index / CHUNK_SIZE
            }
        }
    }

    fn occupancy(&self, chunk: usize) -> i32 {
        self.occupancy[chunk]
    }

    fn is_empty(&self, chunk: usize) -> bool {
        self.occupancy(chunk) == 0
    }

    fn fill_queue(&mut self) {
        self.queue = self
            .counts
            .iter()
            .enumerate()
            .filter_map(|(i, c)| if c == &0 { Some(i) } else { None })
            .collect();
    }

    // Not clear that we really want this kind of logic.
    // This fn may disappear later.
    fn fill_from(&mut self, source: usize, destination: usize) {
        assert_ne!(source, destination);
        // NOTE: code like this could be a perf bottleneck.
        // We may need to get each mutable sub-slice out
        // using .split_at_mut() and some indexing.
        // With the slices, iterators may be helpful
        // in dodging boundary checks.
        for i in 0..self.occupancy(source) {
            self.mutation_ids[destination * CHUNK_SIZE + i as usize] =
                self.mutation_ids[source * CHUNK_SIZE + i as usize];
        }
        self.occupancy[destination] += self.occupancy[source];
    }

    fn position_details<F>(
        &self,
        chunk: usize,
        mutations: &[Mutation],
        f: F,
    ) -> Option<forrustts::Position>
    where
        F: Fn(usize, usize, &[Mutation]) -> forrustts::Position,
    {
        let o = self.occupancy(chunk);
        match o {
            x if x == 0 => None,
            x if x > 0 && ((x as usize) <= CHUNK_SIZE) => Some(f(chunk, x as usize, mutations)),
            _ => panic!("invalid occupancy value"),
        }
    }

    // In general, we don't want to mess w/empty chunks,
    // other than when we FIRST create one.
    // so we may revisit some of this logic later
    // FIXME: code duplication w/last_position
    fn first_position(&self, chunk: usize, mutations: &[Mutation]) -> Option<forrustts::Position> {
        self.position_details(chunk, mutations, |c, _, m| {
            m[self.mutation_ids[c * CHUNK_SIZE] as usize].position()
        })
        //let o = self.occupancy(chunk);
        //match o {
        //    x if x == 0 => None,
        //    x if x > 0 && ((x as usize) <= CHUNK_SIZE) => {
        //        Some(mutations[self.mutation_ids[chunk * CHUNK_SIZE] as usize].position())
        //    }
        //    _ => panic!("invalid occupancy value"),
        //}
    }

    // In general, we don't want to mess w/empty chunks,
    // other than when we FIRST create one.
    // so we may revisit some of this logic later
    // FIXME: code duplication w/first_position
    fn last_position(&self, chunk: usize, mutations: &[Mutation]) -> Option<forrustts::Position> {
        self.position_details(chunk, mutations, |c, x, m| {
            m[self.mutation_ids[c * CHUNK_SIZE + (x - 1)] as usize].position()
        })
        //let o = self.occupancy(chunk);
        //match o {
        //    x if x == 0 => None,
        //    x if x > 0 && ((x as usize) <= CHUNK_SIZE) => Some(
        //        mutations[self.mutation_ids[chunk * CHUNK_SIZE + (x - 1) as usize] as usize]
        //            .position(),
        //    ),
        //    _ => panic!("invalid occupancy value"),
        //}
    }
}

#[derive(Default)]
struct HaploidGenomes {
    mutation_chunk_ids: Vec<u32>, // each genome is a CONTIGUOUS range of chunk indexes
    starts: Vec<usize>,           // For each genome, where is its first chunk?
    stops: Vec<usize>, // One past last chunk such that a genome is mutation_chunk_ids[starts[i]..stops[i]]
}

#[cfg(test)]
mod test_mutation_chunks {
    use super::Mutation;
    use super::MutationChunks;
    use super::CHUNK_SIZE;

    #[test]
    fn test_add_chunk() {
        let mut mc = MutationChunks::default();
        let nc = mc.new_chunk();
        assert_eq!(nc, 0);
        assert_eq!(mc.occupancy(nc), 0);
    }

    #[test]
    fn test_recycle_chunk() {
        let mut mc = MutationChunks::default();
        let nc = mc.new_chunk();
        assert_eq!(nc, 0);
        assert!(mc.is_empty(nc));
        assert!(mc.last_position(nc, &[]).is_none());
        assert!(mc.first_position(nc, &[]).is_none());
        mc.fill_queue();
        let nc = mc.new_chunk();
        assert_eq!(nc, 0);
        mc.counts[0] += 1;
        mc.fill_queue();
        let nc = mc.new_chunk();
        assert_eq!(nc, 1);
        mc.counts[1] += 1;
        mc.fill_queue();
        let nc = mc.new_chunk();
        assert_eq!(nc, 2);
    }

    #[test]
    fn test_first_last_pos_full_chunk() {
        let mut mc = MutationChunks::default();
        let first = mc.new_chunk();

        // Make up some data
        let mut mutations = vec![];
        for i in 0..CHUNK_SIZE {
            mc.mutation_ids[first * CHUNK_SIZE + i] = i.try_into().unwrap();
            let pos = i64::try_from(i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
        }
        mc.occupancy[first] = CHUNK_SIZE as i32;
        assert!(mc.first_position(first, &mutations).is_some());
        assert!(mc.last_position(first, &mutations).is_some());
        assert_eq!(
            mc.first_position(first, &mutations),
            Some(forrustts::Position::try_from(0_i64).unwrap())
        );
        assert_eq!(
            mc.last_position(first, &mutations),
            Some(forrustts::Position::try_from((CHUNK_SIZE as i64) - 1).unwrap())
        );
    }

    #[test]
    fn test_fill_entire_chunk_from_another_full_chunk() {
        let mut mc = MutationChunks::default();
        let first = mc.new_chunk();
        let second = mc.new_chunk();

        // Make up some data
        for i in 0..CHUNK_SIZE {
            mc.mutation_ids[first * CHUNK_SIZE + i] = i.try_into().unwrap();
        }
        mc.occupancy[first] = CHUNK_SIZE as i32;

        mc.fill_from(first, second);
        let s = &mc.mutation_ids[second * CHUNK_SIZE..];
        for i in 0..CHUNK_SIZE {
            assert_eq!(s[i], i as u32);
        }
        assert_eq!(mc.occupancy(second), CHUNK_SIZE as i32);
    }
}
