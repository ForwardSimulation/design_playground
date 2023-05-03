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
    occupancy: Vec<i32>,
}

#[derive(Default)]
struct HaploidGenomes {
    mutation_chunk_ids: Vec<u32>, // each genome is a CONTIGUOUS range of chunk indexes
    starts: Vec<usize>,           // For each genome, where is its first chunk?
    stops: Vec<usize>, // One past last chunk such that a genome is mutation_chunk_ids[starts[i]..stops[i]]
}
