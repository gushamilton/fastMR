// Exact per-exposure greedy LD clumping against an in-memory LD graph.
#include <Rcpp.h>
#include <vector>
#include <cstring>

using namespace Rcpp;

// Split PLINK2 .vcor lines into the ID_A (field 3) and ID_B (field 6) columns
// without allocating a per-line list.  Header and short lines are skipped.
// [[Rcpp::export(name = ".fastmr_vcor_ids")]]
List fastmr_vcor_ids(CharacterVector lines) {
  R_xlen_t n = lines.size();
  std::vector<SEXP> a, b;
  a.reserve(n);
  b.reserve(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (lines[i] == NA_STRING) continue;
    const char* s = CHAR(STRING_ELT(lines, i));
    if (s[0] == '#' || s[0] == '\0') continue;
    const char* tab[6];
    int found = 0;
    const char* p = s;
    const char* start[7];
    start[0] = s;
    for (; *p; ++p) {
      if (*p == '\t') {
        tab[found] = p;
        ++found;
        start[found] = p + 1;
        if (found == 6) break;
      }
    }
    if (found < 5) continue;
    const char* e3 = tab[2];
    // field 3 = (start[2], tab[2]); field 6 = (start[5], tab[5] or end)
    const char* f6e = (found >= 6) ? tab[5] : s + std::strlen(s);
    a.push_back(Rf_mkCharLenCE(start[2], (int)(e3 - start[2]), CE_UTF8));
    b.push_back(Rf_mkCharLenCE(start[5], (int)(f6e - start[5]), CE_UTF8));
  }
  R_xlen_t m = a.size();
  CharacterVector ra(m), rb(m);
  for (R_xlen_t i = 0; i < m; ++i) {
    SET_STRING_ELT(ra, i, a[i]);
    SET_STRING_ELT(rb, i, b[i]);
  }
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
