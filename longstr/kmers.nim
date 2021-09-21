import kmer
import strutils
import math
import tables
import nimpy

{.push checks:off optimization:speed.}
iterator slide_by*(s:string, k: int): uint64 {.inline.} =
  ## given a string (DNA seq) yield the minimum kmer on the forward strand
  if k <= s.len:
    var base: char
    var f = s[0..<k].encode
    var kmin = f
    # note that we are just rotating the kmer here, not adding new bases
    for j in 0..<k:
        base = s[j]
        f.forward_add(base, k)
        kmin = min(f, kmin)
    yield kmin

    # after the first k, then we can use the forward add of k bases
    # to get to the next encode
    for i in countup(k, s.high - k + 1, k):
      for m in 0..<k:
        f.forward_add(s[i + m], k)
      kmin = f
      # then rotate the kmer
      for j in 0..<k:
          base = s[i + j]
          f.forward_add(base, k)
          kmin = min(f, kmin)
      yield kmin
{.pop.}

type Seq*[T] = object
  imax*: int
  A*: seq[T]

type Seqs*[T] = array[7, Seq[T]]

# Report the n most frequent keys in a CountTable, in decending order
proc most_frequent*[T](table: var CountTable[T], n: int): seq[T] =
  table.sort()
  if n > len(table):
    raise newException(IndexError, &"Insufficient keys in CountTable ({len(table)}) to report {n}")
  result = newseq[T](n)
  var i = 0
  for key, count in table:
    if i < n:
      result[i] = key
      inc i
    else:
      break

proc isNaN(v:float): bool {.inline.} =
  return v.classify == fcNaN

proc init*[T](): Seqs[T] =
  result = [
     Seq[T](A: newSeq[T](0)),
     Seq[T](A: newSeq[T](0)),
     Seq[T](A: newSeq[T](16)),
     Seq[T](A: newSeq[T](64)),
     Seq[T](A: newSeq[T](256)),
     Seq[T](A: newSeq[T](1024)),
     Seq[T](A: newSeq[T](4096)),
  ]

proc inc*[T](s:var Seq[T], enc:uint64) {.inline.} =
  s.A[enc.int].inc
  if s.imax == -1 or s.A[enc] > s.A[s.imax]:
    s.imax = enc.int

proc argmax*[T](s: Seq[T]): uint64 {.inline.} =
  return s.imax.uint64

proc clear*[T](s: var Seq[T]) {.inline.} =
  if s.imax == -1: return
  zeroMem(s.A[0].addr, sizeof(T) * len(s.A))
  s.imax = -1

proc count*(read: var string, k: int, count: var Seq[uint8]): int {.inline.} =
  # count the repeats of length k in read and return the most frequent
  count.clear
  for enc in read.slide_by(k):
    count.inc(enc)
  if count.imax == -1: return 0
  result = count.A[count.imax].int

proc reduce_repeat*(rep: var array[6, char]): int =
  ## change e.g. AA to A. result is the mutliplier, so AA to A becomes 2.
  ## CCC to C is 1.
  result = 1
  if rep[0] == '\0': return
  let seen = rep[0]
  for i in 1..<rep.len:
    if rep[i] == '\0': break
    if rep[i] != seen: return
  # if we got here, we have a homopolymer
  for i in 1..<rep.len:
    if rep[i] == '\0': break
    result.inc
    rep[i] = '\0'

proc tostring*(a:array[6, char]): string =
  for c in a:
    if c == 0.char: return
    result.add(c)

proc `<`(a: array[6, char], b: array[6, char]): bool {.inline.} =
  if a[0] != b[0]:
    return a[0] < b[0]
  if a[1] != b[1]:
    return a[1] < b[1]
  if a[2] != b[2]:
    return a[2] < b[2]
  if a[3] != b[3]:
    return a[3] < b[3]
  if a[4] != b[4]:
    return a[4] < b[4]
  return a[5] < b[5]

# This is the bottleneck for run time at the moment
proc top_kmer*(read: var string, motif_count: var int): array[6, char] =
#tuple[motif: string, motif_count: float] =
  let proportion_repeat = 0.0 # Remove this?
  var counts = init[uint8]()
  motif_count = 0
  if read.count('N') > 20: return
  var s = newString(6)

  var best_score: int = -1
  for k in 2..6:
    var count = read.count(k, counts[k])
    s = newString(k)
    counts[k].argmax.decode(s)
    var score = count * k

    if score <= best_score:
      if count < (read.len.float * 0.12 / k.float).int:
        break
      continue
    count = read.count(s)
    score = count * k
    if score < best_score: continue

    best_score = score
    if count > (read.len.float * proportion_repeat / k.float).int:
      # now check the actual string because the kmer method can't track phase
      if count >= (read.len.float * proportion_repeat / k.float).int:
        copyMem(result.addr, s[0].addr, k)
        motif_count = count
        if motif_count > 0 and result[0] == '\0':
          quit "bad:" & $k & " " & $result & " " & "kmer:" & s
    # if we get here, then we could have a count from the kmers of $n for a
    # reduced rotated kmer of, e.g. TG, but the real repat might be
    # GTCGTCGTC... in which case there would be 0 of TG but we want to continue
    # to find GTC

  motif_count *= reduce_repeat(result)

proc kmer_freq*(s: string): tuple[kmer: string, prop: float] {.exportpy.} =
  #echo s.count_kmers()
  var test_seq = s
  var kmer_count = 0
  var top_kmer = test_seq.top_kmer(kmer_count).tostring()
  var top_prop = kmer_count * top_kmer.len() / test_seq.len()
  result = (top_kmer, top_prop)
