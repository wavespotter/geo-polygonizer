# geo-polygonizer

`geo-polygonizer` is a Rust polygonizer built on top of the [`geo`](https://crates.io/crates/geo) ecosystem.

In basic terms, it takes linework (edges/segments) and reconstructs polygon rings from that graph. It is intended for cases where you have boundaries represented as lines and want valid polygon geometry as output.

## Background

This project is inspired by polygonizer implementations in tools like JTS, GEOS, and Turf, but it is **not** a literal port of any one implementation.

The algorithms and behavior are adapted for this Rust codebase and its dependencies.

## Development

Run tests:

```bash
cargo test
```

Check formatting:

```bash
cargo fmt --all -- --check
```

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
