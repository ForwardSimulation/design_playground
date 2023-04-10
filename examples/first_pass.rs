use clap::Parser;
use forrustts::genetics::GeneticMap;
use forrustts::genetics::GeneticMapBuilder;
use forrustts::genetics::PoissonCrossover;

use design_playground::first_pass::*;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long)]
    seed: u64,
    #[arg(short, long)]
    popsize: u32,
    #[arg(short, long)]
    mu: f64,
    #[arg(short, long)]
    r: f64,
    #[arg(short, long)]
    nsteps: u32,
}

fn main() {
    let args = Args::parse();
    let params = SimParams {
        seed: args.seed,
        size: args.popsize,
        num_generations: args.nsteps,
        mutation_rate: args.mu,
        recrate: args.r,
    };
    let poisson = vec![PoissonCrossover::new(0, 1000000, params.recrate).unwrap()];
    let builder = GeneticMapBuilder::default().extend_poisson(&poisson);
    let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
    evolve_pop_with_haplotypes(params, genetic_map).unwrap();
}
