# Finite-coupling spectrum of O(N) model in AdS — `ONinAdS` package

Accompanying code for "Finite-coupling spectrum of O(N) model in AdS" [2503.16345 \[hep-th\]](https://arxiv.org/abs/2503.16345).

## Software Requirements

- [Wolfram Mathematica](https://www.wolfram.com/mathematica/) — preferably version `>=14.0`

## How to setup `ONinAdS` package

We need to copy/symlink the package directory to place where `Mathematica` searches
for packages. This is any directory in Mathematica `$PATH`, for example something
like `/home/username/Wolfram/Applications`.
On `linux`, this can be done using
```sh
ln -sv -t /path/to/directory/in/Mathematica/$Path /path/to/repository-ONinAdS/ONinAdS
```

Now we can load the whole package with (note the trailing backtick)
```mma
Get["ONinAdS`"]
(* or alternatively with << ONinAdS` *)
```

Alternatively, one can provide the full path to the package as
```mma
Get["/path/to/repository-ONinAdS/ONinAdS`"]
```

## Package Showcase

See the `ONinAdS-showcase.nb` file for a quick example of how to use the package.
However, most of the interesting details are directly included in the package files (located in the `ONinAdS/` directory).

## Presentations

In this repository we also include the slides from the following presentations:
- March 27th 2025 — [ONinAdS-talk.pdf](https://jdujava.github.io/ONinAdS/Presentations/ONinAdS-talk.pdf)
- April 4th 2025 — [ONinAdS-seminar.pdf](https://jdujava.github.io/ONinAdS/Presentations/ONinAdS-seminar.pdf)
