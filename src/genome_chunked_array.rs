// Mutations will be stored in blocks of 8
// indexes
static CHUNK_SIZE: u32 = 8;

struct Chunk {
    start: u32,
    stop: u32,
}

struct MutationChunks {
    mutations: Vec<u32>,
    chunks: Vec<rclite::Rc<Chunk>>,
}

struct Individuals {
    chunks: Vec<rclite::Rc<Chunk>>,
    starts: Vec<usize>,
    stops: Vec<usize>,
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
