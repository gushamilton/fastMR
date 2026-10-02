// Exact per-exposure greedy LD clumping against an in-memory LD graph.
#include <Rcpp.h>
#include <vector>
#include <cstring>
#include <cstdio>
#include <memory>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

// Parse an uncompressed PLINK2 .vcor table straight into 1-based vertex ids.
// `ids` are the candidate SNP IDs (vertex v = position v in `ids`).  The first
// non-empty line must be a '#' header naming ID_A before ID_B; any other first
// line is an error.  Blank lines and later '#' lines are skipped; a final line
// without a newline is kept; a trailing '\r' is ignored.  A data line with too
// few fields or an ID outside `ids` is an error.  A zero-byte file is an empty
// graph.  Only std containers are used until the result vectors are built, so
// there is nothing for the GC to trip over.
// [[Rcpp::export(name = ".fastmr_vcor_read")]]
List fastmr_vcor_read(std::string path, CharacterVector ids) {
  std::unordered_map<std::string, int> index;
  index.reserve((size_t)ids.size() * 2 + 16);
  for (R_xlen_t i = 0; i < ids.size(); ++i) {
    if (ids[i] == NA_STRING) stop("NA candidate SNP id");
    index.emplace(std::string(CHAR(STRING_ELT(ids, i))), (int)(i + 1));
  }
  std::unique_ptr<FILE, int (*)(FILE*)> fp(std::fopen(path.c_str(), "rb"), &std::fclose);
  if (!fp) stop("cannot open PLINK2 .vcor file: " + path);
  std::vector<int> lead, target;
  const size_t CH = 1u << 22;
  std::vector<char> buf(CH + 1);
  size_t have = 0;  // carried bytes of a partial line at buf[0..have)
  int fa = 0, fb = 0;
  bool header_done = false;
  long long short_lines = 0;
  std::string key;
  std::string bad;
  auto handle = [&](const char* s, const char* e) {
    if (e > s && e[-1] == '\r') --e;
    if (s == e) return;
    if (!header_done) {
      header_done = true;
      if (*s == '#') {
        int ia = 0, ib = 0, k = 0;
        const char* f = s + 1;
        while (true) {
          const char* t = f;
          while (t < e && *t != '\t') ++t;
          ++k;
          size_t len = (size_t)(t - f);
          if (!ia && len == 4 && !std::strncmp(f, "ID_A", 4)) ia = k;
          if (!ib && len == 4 && !std::strncmp(f, "ID_B", 4)) ib = k;
          if (t >= e) break;
          f = t + 1;
        }
        if (ia && ib && ia < ib) { fa = ia; fb = ib; }
      }
      if (!fa) stop("unrecognised PLINK2 .vcor header (need ID_A and ID_B): " +
                    std::string(s, std::min<size_t>((size_t)(e - s), 200)));
      return;
    }
    if (*s == '#') return;
    const char* sa = nullptr; const char* ea = nullptr;
    const char* sb = nullptr; const char* eb = nullptr;
    int k = 1;
    const char* f = s;
    while (true) {
      const char* t = f;
      while (t < e && *t != '\t') ++t;
      if (k == fa) { sa = f; ea = t; }
      if (k == fb) { sb = f; eb = t; break; }
      if (t >= e) break;
      f = t + 1;
      ++k;
    }
    if (!sb) { ++short_lines; return; }
    key.assign(sa, ea);
    auto ita = index.find(key);
    if (ita == index.end()) stop("PLINK2 .vcor ID not among the candidate SNPs: " + key);
    key.assign(sb, eb);
    auto itb = index.find(key);
    if (itb == index.end()) stop("PLINK2 .vcor ID not among the candidate SNPs: " + key);
    lead.push_back(ita->second);
    target.push_back(itb->second);
  };
  while (true) {
    size_t got = std::fread(buf.data() + have, 1, CH - have, fp.get());
    if (got == 0) {
      if (std::ferror(fp.get())) stop("error reading PLINK2 .vcor file");
      break;
    }
    size_t n = have + got;
    size_t pos = 0;
    while (true) {
      const char* nl = (const char*)std::memchr(buf.data() + pos, '\n', n - pos);
      if (!nl) break;
      handle(buf.data() + pos, nl);
      pos = (size_t)(nl - buf.data()) + 1;
    }
    have = n - pos;
    if (have == CH) stop("PLINK2 .vcor line longer than the read buffer");
    if (have) std::memmove(buf.data(), buf.data() + pos, have);
  }
  if (have) handle(buf.data(), buf.data() + have);  // unterminated final line
  if (short_lines)
    stop("malformed PLINK2 .vcor line(s): " + std::to_string(short_lines) +
         " line(s) have fewer than " + std::to_string(fb) + " fields");
  R_xlen_t m = (R_xlen_t)lead.size();
  IntegerVector ra(m), rb(m);
  std::copy(lead.begin(), lead.end(), ra.begin());
  std::copy(target.begin(), target.end(), rb.begin());
  return List::create(_["lead"] = ra, _["target"] = rb);
}

// Greedy clump of every exposure on one symmetric LD graph.
//   n_snp      number of graph vertices
//   ea, eb     0-based edge endpoints (each undirected edge once or twice)
//   row_snp    0-based vertex of each candidate row, grouped by exposure and
//              ALREADY in greedy (p, SNP) order within each exposure
//   exp_start  0-based offsets into row_snp, length n_exposure + 1
// Within an exposure vertices must be unique.  Returns TRUE for retained rows.
// [[Rcpp::export(name = ".fastmr_graph_clump")]]
LogicalVector fastmr_graph_clump(int n_snp, IntegerVector ea, IntegerVector eb,
                                 IntegerVector row_snp, IntegerVector exp_start) {
  R_xlen_t ne = ea.size();
  if (eb.size() != ne) stop("edge vectors differ in length");
  std::vector<int> deg(n_snp + 1, 0);
  for (R_xlen_t i = 0; i < ne; ++i) {
    int a = ea[i], b = eb[i];
    if (a < 0 || b < 0 || a >= n_snp || b >= n_snp) stop("edge endpoint out of range");
    if (a == b) continue;
    ++deg[a + 1];
    ++deg[b + 1];
  }
  for (int v = 0; v < n_snp; ++v) deg[v + 1] += deg[v];
  std::vector<int> adj(deg[n_snp]);
  std::vector<int> fill(deg.begin(), deg.end() - 1);
  for (R_xlen_t i = 0; i < ne; ++i) {
    int a = ea[i], b = eb[i];
    if (a == b) continue;
    adj[fill[a]++] = b;
    adj[fill[b]++] = a;
  }
  R_xlen_t nrow = row_snp.size();
  LogicalVector keep(nrow);
  std::vector<int> stamp(n_snp, -1), loc(n_snp, 0);
  int n_exp = exp_start.size() - 1;
  std::vector<char> dead;
  for (int e = 0; e < n_exp; ++e) {
    int s = exp_start[e], t = exp_start[e + 1];
    if (s < 0 || t > nrow || t < s) stop("invalid exposure offsets");
    dead.assign(t - s, 0);
    for (int k = s; k < t; ++k) {
      int v = row_snp[k];
      if (v < 0 || v >= n_snp) stop("row vertex out of range");
      stamp[v] = e;
      loc[v] = k - s;
    }
    for (int k = s; k < t; ++k) {
      if (dead[k - s]) continue;
      keep[k] = true;
      int v = row_snp[k];
      for (int j = deg[v]; j < deg[v + 1]; ++j) {
        int w = adj[j];
        if (stamp[w] == e && loc[w] > k - s) dead[loc[w]] = 1;
      }
    }
  }
  return keep;
}
