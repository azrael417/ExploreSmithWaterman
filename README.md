# Reads & References

Load a read and try to align it with different reference sequences.

1. The initial use case we target is 1 read to many references
1. sequences could be 50 to 100,000 characters long
1. all pairs of (read, reference) are independent

```c++
  while (fastq.LoadNextRead(&readname, &sequence, &qual)) {
    for (int i = 0; i < refs_count; ++i) {
      const char* pReference = refs.GetReferenceSequence(i, &length);
      sw.Align(&alignment, cigarSW, pReference, length, sequence.c_str(), sequence.size());
      PrintAlignment(readname, sequence, cigarSW, alignment);
    }
  }
```
  
# Alignment
  
Alignment is done with [Smith-Waterman](https://en.wikipedia.org/wiki/Smith–Waterman_algorithm).
  
1. Compute a scoring matrix
1. Trackback to find the best alignment

- full of conditionals
- dependencies in scoring matrix `H(i,j)` depends on `H(i-1,j-1)`, `H(i-1,j)` & `H(i,j-1)`

# References

- [SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in Genomic Applications](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138)

# Changelog 

1. pull memory (re)allocation's from inner loop
1. option to load new reads in batches
1. turned everything into views
