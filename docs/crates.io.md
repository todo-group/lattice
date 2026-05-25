# Publishing to crates.io

The primary crates.io package is `lattice-core`.

`lattice-python` is a PyPI extension crate and should normally not be published
to crates.io. `lattice-ffi` is optional; publish it only if downstream Rust/C
consumers need the C ABI crate directly.

## One-time setup

Create a crates.io account, verify your email address, and create an API token:

```text
https://crates.io/me
```

Store the token locally:

```sh
cargo login
```

The token is secret. Revoke it on crates.io if it is ever exposed.

## Preflight checks

Run tests:

```sh
cargo test -p lattice-core
```

Check the package without uploading:

```sh
cargo publish -p lattice-core --dry-run
```

Inspect the packaged files:

```sh
cargo package -p lattice-core --list
```

If the package contains unnecessary files, add `include` or `exclude` entries to
`rust/lattice-core/Cargo.toml` before publishing.

## Publish lattice-core

Publish:

```sh
cargo publish -p lattice-core
```

After publishing, docs.rs should build the documentation automatically:

```text
https://docs.rs/lattice-core
```

## Publishing lattice-ffi, if needed

Publish `lattice-ffi` only after `lattice-core` is available on crates.io.

Before publishing, make sure the dependency has both a path for local workspace
development and a version for crates.io resolution:

```toml
lattice-core = { version = "0.1.0", path = "../lattice-core" }
```

Then run:

```sh
cargo test -p lattice-ffi
cargo publish -p lattice-ffi --dry-run
cargo package -p lattice-ffi --list
cargo publish -p lattice-ffi
```

## Versioning rules

crates.io releases are permanent:

* a published version cannot be overwritten
* uploaded package contents cannot be deleted
* fixes require a new version

For example, after publishing `0.1.0`, a fix should be released as `0.1.1` or a
later SemVer-compatible version.

## Release checklist

1. Update the crate version in `Cargo.toml`.
2. Run `cargo test -p lattice-core`.
3. Run `cargo publish -p lattice-core --dry-run`.
4. Inspect `cargo package -p lattice-core --list`.
5. Commit the release changes.
6. Tag the release, for example `v0.1.0`.
7. Run `cargo publish -p lattice-core`.
8. Confirm the crate page and docs.rs page are live.
