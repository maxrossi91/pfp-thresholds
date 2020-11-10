# pfp-thresholds
Computing thresholds from Prefix-free parsing

## Download

```console
git clone https://github.com/maxrossi91/pfp-thresholds
```

## Compile

```console
mkdir build
cd build; cmake ..
make
```

### Debug

```console
mkdir build
cd build; cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

## Contents

### LCP computation

* `pfp_lcp`: build two arrays reporting, for each run [l..r] of the BWT, the position and the value of the minimum LCP[l..r+1] and its position. The two arrays are built from the prefix-free parsing.

### Thresholds computation

* `pfp_thrersolds`: build the thresholds from the prefix-free parsing. (`BigBWT`, `pfp_thresholds.cpp`)

* `gsacak_thresholds`: build the thresholds using `gsacak`. (`gsacak`)

* `bwt_lcp_thresholds`: build the thresholds using `RLBWT2LCP`. (`BigBWT`,`RLBWT2LCP`)

### Matching Statistics

* `matching_statistics`: computes the matching statistics from the BWT and the thresholds, using the parsing for random access.

* `sdsl_matching_statistics`: computes the matching statistics from the text using `sdsl`.

## Authors 

### Theoretical results:

* Christina Boucher
* Travis Gagie
* Massimiliano Rossi

### Implementation:

* Massimiliano Rossi

### Experiments:

* Massimiliano Rossi
* Marco Oliva