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

    fn position_details<F>(&self, chunk: usize, f: F) -> Option<forrustts::Position>
    where
        F: Fn(usize) -> forrustts::Position,
    {
        let o = self.occupancy(chunk);
        match o {
            x if x == 0 => None,
            x if x > 0 && ((x as usize) <= CHUNK_SIZE) => Some(f(x as usize)),
            _ => panic!("invalid occupancy value"),
        }
    }

    // In general, we don't want to mess w/empty chunks,
    // other than when we FIRST create one.
    // so we may revisit some of this logic later
    fn first_position(&self, chunk: usize, mutations: &[Mutation]) -> Option<forrustts::Position> {
        self.position_details(chunk, |_| {
            mutations[self.mutation_ids[chunk * CHUNK_SIZE] as usize].position()
        })
    }

    // In general, we don't want to mess w/empty chunks,
    // other than when we FIRST create one.
    // so we may revisit some of this logic later
    fn last_position(&self, chunk: usize, mutations: &[Mutation]) -> Option<forrustts::Position> {
        self.position_details(chunk, |x| {
            mutations[self.mutation_ids[chunk * CHUNK_SIZE + (x - 1)] as usize].position()
        })
    }
}

// This is a quick implementation of C++'s remove_if, which is a method
// of partitioning a slice such that retained items are STABLY moved to the
// front.
// NOTE: we can probably do this with generics later on
// NOTE: should the be associated w/HaploidGenomes?
fn remove_chunk_if_empty(
    mutation_chunk_ids: &mut [u32],
    mutation_chunks: &MutationChunks,
) -> usize {
    let mut first = match mutation_chunk_ids
        .iter()
        .position(|&c| mutation_chunks.is_empty(c as usize))
    {
        Some(index) => index,
        None => mutation_chunk_ids.len(),
    };
    let f = first;
    for i in f..mutation_chunk_ids.len() {
        if !mutation_chunks.is_empty(mutation_chunk_ids[i] as usize) {
            mutation_chunk_ids.swap(first, i);
            //let temp = mutation_chunk_ids[first];
            //mutation_chunk_ids[first] = mutation_chunk_ids[i];
            //mutation_chunk_ids[i] = temp;
            first += 1;
        }
    }
    first
}

#[derive(Default)]
struct HaploidGenomes {
    mutation_chunk_ids: Vec<u32>, // each genome is a CONTIGUOUS range of chunk indexes
    starts: Vec<usize>,           // For each genome, where is its first chunk?
    stops: Vec<usize>, // One past last chunk such that a genome is mutation_chunk_ids[starts[i]..stops[i]]
}

impl HaploidGenomes {
    fn move_empty_chunks(&mut self, mutation_chunks: &MutationChunks) {
        assert_eq!(self.starts.len(), self.stops.len());
        let ngenomes = self.starts.len();
        for genome in 0..ngenomes {
            let chunks = &mut self.mutation_chunk_ids[self.starts[genome]..self.stops[genome]];
            let nremoved = remove_chunk_if_empty(chunks, mutation_chunks);
            self.stops[0] = nremoved;
        }
    }
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

#[cfg(test)]
mod test_haploid_genomes {
    use super::*;

    // I STRONGLY DISLIKE the semantics here.
    // This will lead to all sorts of edge cases
    #[test]
    fn test_search_thru_empty_chunk() {
        let mut mc = MutationChunks::default();
        let first = mc.new_chunk();
        let second = mc.new_chunk();
        let third = mc.new_chunk();
        let fourth = mc.new_chunk();
        let mut mutations = vec![];
        for i in 0..CHUNK_SIZE {
            mc.mutation_ids[first * CHUNK_SIZE + i] = i.try_into().unwrap();
            let pos = i64::try_from(i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
        }
        for i in 0..CHUNK_SIZE {
            let pos = i64::try_from(CHUNK_SIZE + i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
            mc.mutation_ids[third * CHUNK_SIZE + i] = (mutations.len() - 1).try_into().unwrap();
        }
        mc.occupancy[first] = CHUNK_SIZE as i32;
        mc.occupancy[third] = CHUNK_SIZE as i32;

        assert!(mc.first_position(second, &mutations).is_none());
        assert!(mc.last_position(second, &mutations).is_none());
        assert!(mc.first_position(fourth, &mutations).is_none());
        assert!(mc.last_position(fourth, &mutations).is_none());
        assert_eq!(
            mc.first_position(first, &mutations),
            Some(forrustts::Position::try_from(0).unwrap())
        );
        assert_eq!(
            mc.last_position(first, &mutations),
            Some(forrustts::Position::try_from((CHUNK_SIZE as i64) - 1).unwrap())
        );
        assert_eq!(
            mc.first_position(third, &mutations),
            Some(forrustts::Position::try_from(CHUNK_SIZE as i64).unwrap())
        );
        assert_eq!(
            mc.last_position(third, &mutations),
            Some(forrustts::Position::try_from((2 * CHUNK_SIZE as i64) - 1).unwrap())
        );

        // Okay, so now we can set up some haploid genomes
        let mut haploid_genomes = HaploidGenomes::default();
        haploid_genomes.mutation_chunk_ids.push(first as u32);
        haploid_genomes.mutation_chunk_ids.push(second as u32);
        haploid_genomes.mutation_chunk_ids.push(third as u32);
        haploid_genomes.mutation_chunk_ids.push(fourth as u32);
        haploid_genomes.starts.push(0);
        haploid_genomes.stops.push(4);

        // Okay, now we have a genome with chunks that are full, empty, full.
        // We envision this happening when, at the end of a generation,
        // all mutations in chunk 1 (2nd chunk) become fixed. The idea is that
        // we will mark all elements corresponding to fixed mutations with
        // our sentinel value.

        // The goal of this test is to be able to correctly use partition_point
        // for this kind of case.
        let test_position = forrustts::Position::try_from(100).unwrap();
        assert!(test_position >= mc.first_position(third, &mutations).unwrap());
        assert!(test_position < mc.last_position(third, &mutations).unwrap());

        let genome = &haploid_genomes.mutation_chunk_ids
            [haploid_genomes.starts[0]..haploid_genomes.stops[0]];
        let p = genome.partition_point(|&c| {
            let comp = mc.last_position(c as usize, &mutations).is_none()
                || mc.last_position(c as usize, &mutations).unwrap() < test_position;
            comp
        });
        assert!(p > 0);
        assert_eq!(p, 2);

        let test_position = forrustts::Position::try_from(200).unwrap();
        let p = genome.partition_point(|&c| {
            let comp = mc.last_position(c as usize, &mutations).is_none()
                || mc.last_position(c as usize, &mutations).unwrap() < test_position;
            comp
        });
        assert_eq!(p, 4);
    }

    // This leads to saner/safer semantics than the above
    #[test]
    fn purge_empty_chunks() {
        let mut mc = MutationChunks::default();
        let first = mc.new_chunk();
        let second = mc.new_chunk();
        let third = mc.new_chunk();
        let fourth = mc.new_chunk();
        let fifth = mc.new_chunk();
        let mut mutations = vec![];
        for i in 0..CHUNK_SIZE {
            mc.mutation_ids[first * CHUNK_SIZE + i] = i.try_into().unwrap();
            let pos = i64::try_from(i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
        }
        for i in 0..CHUNK_SIZE {
            let pos = i64::try_from(CHUNK_SIZE + i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
            mc.mutation_ids[third * CHUNK_SIZE + i] = (mutations.len() - 1).try_into().unwrap();
        }
        for i in 0..CHUNK_SIZE {
            let pos = i64::try_from(2 * CHUNK_SIZE + i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
            mc.mutation_ids[fifth * CHUNK_SIZE + i] = (mutations.len() - 1).try_into().unwrap();
        }
        mc.occupancy[first] = CHUNK_SIZE as i32;
        mc.occupancy[third] = CHUNK_SIZE as i32;
        mc.occupancy[fifth] = CHUNK_SIZE as i32;

        assert!(mc.first_position(second, &mutations).is_none());
        assert!(mc.last_position(second, &mutations).is_none());
        assert!(mc.first_position(fourth, &mutations).is_none());
        assert!(mc.last_position(fourth, &mutations).is_none());
        assert_eq!(
            mc.first_position(first, &mutations),
            Some(forrustts::Position::try_from(0).unwrap())
        );
        assert_eq!(
            mc.last_position(first, &mutations),
            Some(forrustts::Position::try_from((CHUNK_SIZE as i64) - 1).unwrap())
        );
        assert_eq!(
            mc.first_position(third, &mutations),
            Some(forrustts::Position::try_from(CHUNK_SIZE as i64).unwrap())
        );
        assert_eq!(
            mc.last_position(third, &mutations),
            Some(forrustts::Position::try_from((2 * CHUNK_SIZE as i64) - 1).unwrap())
        );
        assert_eq!(
            mc.first_position(fifth, &mutations),
            Some(forrustts::Position::try_from(2 * CHUNK_SIZE as i64).unwrap())
        );
        assert_eq!(
            mc.last_position(fifth, &mutations),
            Some(forrustts::Position::try_from((3 * CHUNK_SIZE as i64) - 1).unwrap())
        );

        // Okay, so now we can set up some haploid genomes
        let mut haploid_genomes = HaploidGenomes::default();
        haploid_genomes.mutation_chunk_ids.push(first as u32);
        haploid_genomes.mutation_chunk_ids.push(second as u32);
        haploid_genomes.mutation_chunk_ids.push(third as u32);
        haploid_genomes.mutation_chunk_ids.push(fourth as u32);
        haploid_genomes.mutation_chunk_ids.push(fifth as u32);
        haploid_genomes.starts.push(0);
        haploid_genomes.stops.push(5);

        haploid_genomes.move_empty_chunks(&mc);
        assert_eq!(haploid_genomes.stops[0], 3);

        let genome = &haploid_genomes.mutation_chunk_ids
            [haploid_genomes.starts[0]..haploid_genomes.stops[0]];
        assert_eq!(genome, &[first as u32, third as u32, fifth as u32]);
        assert_eq!(genome.len(), 3);
        assert!(genome
            .windows(2)
            .all(|w| mutations[w[0] as usize].position() <= mutations[w[1] as usize].position()));
        let test_position = forrustts::Position::try_from(100).unwrap();
        let p = genome.partition_point(|&c| {
            let comp = mc.last_position(c as usize, &mutations).unwrap() < test_position;
            comp
        });
        assert!(p > 0);
        assert!(p == 1);

        let test_position = forrustts::Position::try_from(200).unwrap();
        let p = genome.partition_point(|&c| {
            let comp = mc.last_position(c as usize, &mutations).unwrap() < test_position;
            comp
        });
        assert!(p > 0);
        assert!(p == 3);
    }
}

#[cfg(test)]
mod tdd_crossover_semantics {
    use forrustts::Position;

    use super::*;

    fn single_crossover(
        genomes: (usize, usize),
        breakpoint: Position,
        mutations: &[Mutation],
        mutation_chunks: &MutationChunks,
        haploid_genomes: &HaploidGenomes,
        output: &mut Vec<u32>,
    ) {
        let genome0 = &haploid_genomes.mutation_chunk_ids
            [haploid_genomes.starts[genomes.0]..haploid_genomes.stops[genomes.0]];
        let genome1 = &haploid_genomes.mutation_chunk_ids
            [haploid_genomes.starts[genomes.1]..haploid_genomes.stops[genomes.1]];

        let p = genome0.partition_point(|&chunk| {
            let comp = mutation_chunks
                .last_position(chunk as usize, &mutations)
                .unwrap()
                < breakpoint;
            comp
        });
        output.extend_from_slice(&genome0[0..p]);
        let p = genome1.partition_point(|&chunk| {
            let comp = mutation_chunks
                .last_position(chunk as usize, &mutations)
                .unwrap()
                < breakpoint;
            comp
        });
        output.extend_from_slice(&genome1[p..]);
    }

    #[test]
    fn test_simple_merge() {
        let mut mutation_chunks = MutationChunks::default();
        let first = mutation_chunks.new_chunk();
        let second = mutation_chunks.new_chunk();
        let third = mutation_chunks.new_chunk();
        let mut mutations = vec![];
        for i in 0..CHUNK_SIZE {
            mutation_chunks.mutation_ids[first * CHUNK_SIZE + i] = i.try_into().unwrap();
            let pos = i64::try_from(i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
        }
        for i in 0..CHUNK_SIZE {
            let pos = i64::try_from(CHUNK_SIZE + i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
            mutation_chunks.mutation_ids[second * CHUNK_SIZE + i] =
                (mutations.len() - 1).try_into().unwrap();
        }
        for i in 0..CHUNK_SIZE {
            let pos = i64::try_from(2 * CHUNK_SIZE + i).unwrap();
            mutations.push(Mutation::new(pos.try_into().unwrap(), vec![], 0.into()));
            mutation_chunks.mutation_ids[third * CHUNK_SIZE + i] =
                (mutations.len() - 1).try_into().unwrap();
        }
        mutation_chunks.occupancy.fill(CHUNK_SIZE as i32);
        assert_eq!(mutations.len(), 3 * CHUNK_SIZE);
        let mut haploid_genomes = HaploidGenomes::default();

        // make first genome
        haploid_genomes.mutation_chunk_ids.push(first as u32);
        haploid_genomes.mutation_chunk_ids.push(second as u32);
        haploid_genomes.starts.push(0);
        haploid_genomes.stops.push(2);
        // make second genome
        haploid_genomes.mutation_chunk_ids.push(first as u32);
        haploid_genomes.mutation_chunk_ids.push(third as u32);
        haploid_genomes.starts.push(2);
        haploid_genomes.stops.push(4);

        let breakpoint = Position::try_from(CHUNK_SIZE as i64).unwrap();
        let mut output = vec![];
        single_crossover(
            (0, 1),
            breakpoint,
            &mutations,
            &mutation_chunks,
            &haploid_genomes,
            &mut output,
        );
        println!("{output:?}");
        assert_eq!(output, &[0, 2]);

        let breakpoint = Position::try_from(2 * CHUNK_SIZE as i64).unwrap();
        let mut output = vec![];
        single_crossover(
            (0, 1),
            breakpoint,
            &mutations,
            &mutation_chunks,
            &haploid_genomes,
            &mut output,
        );
        println!("{output:?}");
        assert_eq!(output, &[0, 1]);

        let breakpoint = Position::try_from(10 + CHUNK_SIZE as i64).unwrap();
        let mut output = vec![];
        single_crossover(
            (0, 1),
            breakpoint,
            &mutations,
            &mutation_chunks,
            &haploid_genomes,
            &mut output,
        );
        println!("{output:?}");
        assert_eq!(output, &[0, 3]);
    }
}
