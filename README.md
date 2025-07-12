# Finite-coupling spectrum of O(N) model in AdS — Mathematica package `ONinAdS`

The main purpose of this repository is to provide the accompanying code for the preprint:<br>
  Jonáš Dujava and Petr Vaško, [Finite-coupling spectrum of O(N) model in AdS](https://arxiv.org/pdf/2503.16345)<br>
    `arXiv`: [2503.16345 \[hep-th\]](https://arxiv.org/abs/2503.16345), `GitHub`: [jdujava/ONinAdS](https://github.com/jdujava/ONinAdS)

In addition, you can also find here the revised version of my master thesis:<br>
  Jonáš Dujava, [Strongly Coupled Quantum Field Theory in Anti-de Sitter Spacetime](https://jdujava.github.io/ONinAdS/SCQFTinAdS.pdf) — Master Thesis<br>
    `arXiv`: [2507.07111 \[hep-th\]](https://arxiv.org/abs/2507.07111), `GitHub`: [jdujava/ONinAdS](https://github.com/jdujava/ONinAdS)

### Quick Preview
<img width="2776" height="2105" alt="thesis-dual00" src="https://github.com/user-attachments/assets/5c7608e8-4ed1-4e6b-b201-335fbe7fe7e3" />
<img width="2776" height="2105" alt="thesis-dual01" src="https://github.com/user-attachments/assets/376a5b62-5b14-483c-87ba-da0808c51a32" />
<img width="2776" height="2105" alt="thesis-dual03" src="https://github.com/user-attachments/assets/bf719391-f5e7-4e56-baa8-05fe330bbc86" />
<img width="2776" height="2105" alt="thesis-dual04" src="https://github.com/user-attachments/assets/b9ed8eaa-23b9-4cda-87b8-4dcc2510910b" />
<img width="2776" height="2105" alt="thesis-dual05" src="https://github.com/user-attachments/assets/50cd1894-b374-4206-9d81-78498f15c25b" />
<img width="2776" height="2105" alt="thesis-dual06" src="https://github.com/user-attachments/assets/702ce869-1f9c-4ceb-a1ce-f7ce815bfb3f" />
<img width="2776" height="2105" alt="thesis-dual07" src="https://github.com/user-attachments/assets/502eafd8-7c24-47f8-84a2-3fef5415e351" />
<img width="2776" height="2105" alt="thesis-dual08" src="https://github.com/user-attachments/assets/0bd1db62-d2a4-4ff0-b126-0c176ac76123" />
<img width="2776" height="2105" alt="thesis-dual09" src="https://github.com/user-attachments/assets/a10b0779-4fb0-444b-b360-99545077d92a" />

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

This repository also includes the slides from the following presentations:
- March 27th 2025 — [ONinAdS-talk.pdf](https://jdujava.github.io/ONinAdS/Presentations/ONinAdS-talk.pdf)
- April 4th 2025 — [ONinAdS-seminar.pdf](https://jdujava.github.io/ONinAdS/Presentations/ONinAdS-seminar.pdf)
