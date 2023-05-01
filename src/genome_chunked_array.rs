// Mutations will be stored in blocks of 64
// indexes
static CHUNK_SIZE: u32 = 64;

struct Chunk {
    start: u32,
    stop: u32,
}

struct MutationChunks {
    mutations: Vec<u32>,
    chunks: Vec<Chunk>,
}

struct Genomes {
    chunks: Vec<u32>,   // each genome is a range of chunk indexes
    starts: Vec<usize>, // For each genome, where is its first chunk?
    stops: Vec<usize>,  // One past last chunk such that a genome is chunks[starts[i]..stops[i]]
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
}
