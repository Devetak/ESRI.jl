#include <armadillo>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

arma::sp_mat scale_each_row(const arma::sp_mat& input, const arma::vec& factors) {
    arma::sp_mat output(input);
    arma::sp_mat::const_iterator start = input.begin();
    arma::sp_mat::const_iterator end = input.end();
    for (arma::sp_mat::const_iterator it = start; it != end; ++it) {
        output(it.row(), it.col()) *= factors(it.row());
    }
    return output;
}

arma::uvec read_uvec_1based(const std::string& path) {
    std::ifstream io(path);
    if (!io) {
        throw std::runtime_error("failed to open " + path);
    }
    std::vector<arma::uword> values;
    arma::uword value = 0;
    while (io >> value) {
        values.push_back(value - 1);
    }
    arma::uvec out(values.size());
    for (arma::uword i = 0; i < values.size(); ++i) {
        out(i) = values[i];
    }
    return out;
}

arma::uvec read_uvec_raw(const std::string& path) {
    std::ifstream io(path);
    if (!io) {
        throw std::runtime_error("failed to open " + path);
    }
    std::vector<arma::uword> values;
    arma::uword value = 0;
    while (io >> value) {
        values.push_back(value);
    }
    arma::uvec out(values.size());
    for (arma::uword i = 0; i < values.size(); ++i) {
        out(i) = values[i];
    }
    return out;
}

arma::sp_mat read_edges(const std::string& path, arma::uword n) {
    std::ifstream io(path);
    if (!io) {
        throw std::runtime_error("failed to open " + path);
    }

    std::vector<arma::uword> rows;
    std::vector<arma::uword> cols;
    std::vector<double> vals;
    arma::uword i = 0;
    arma::uword j = 0;
    double x = 0.0;
    while (io >> i >> j >> x) {
        rows.push_back(i - 1);
        cols.push_back(j - 1);
        vals.push_back(x);
    }

    arma::umat locations(2, rows.size());
    arma::vec values(vals.size());
    for (arma::uword k = 0; k < rows.size(); ++k) {
        locations(0, k) = rows[k];
        locations(1, k) = cols[k];
        values(k) = vals[k];
    }
    return arma::sp_mat(locations, values, n, n);
}

arma::sp_mat build_psup(const arma::uvec& p_cons, arma::uword m) {
    arma::umat locations(2, p_cons.n_elem);
    arma::vec values(p_cons.n_elem, arma::fill::ones);
    for (arma::uword i = 0; i < p_cons.n_elem; ++i) {
        locations(0, i) = i;
        locations(1, i) = p_cons(i);
    }
    return arma::sp_mat(locations, values, p_cons.n_elem, m);
}

arma::sp_mat row_normalize(const arma::sp_mat& input, const arma::vec& denom) {
    arma::sp_mat output(input);
    for (arma::sp_mat::iterator it = output.begin(); it != output.end(); ++it) {
        const double d = denom(it.row());
        *it = (d == 0.0) ? 0.0 : (*it / d);
    }
    return output;
}

arma::sp_mat lambda_d_calc(
    const arma::sp_mat& W,
    const arma::uvec& p_cons,
    const arma::uvec& essential_flags
) {
    const arma::uword n = W.n_rows;
    const arma::uword m = essential_flags.n_elem;
    arma::sp_mat psup = build_psup(p_cons, m);
    arma::mat pi_abs = arma::mat(arma::trans(W) * psup);
    arma::vec in_str = arma::vec(arma::mat(arma::sum(W, 0).t()));

    arma::sp_mat lambda(W);
    for (arma::sp_mat::iterator it = lambda.begin(); it != lambda.end(); ++it) {
        const arma::uword supplier = it.row();
        const arma::uword customer = it.col();
        const arma::uword sector = p_cons(supplier);
        double denom = 0.0;
        if (essential_flags(sector) == 1) {
            denom = pi_abs(customer, sector);
        } else {
            denom = in_str(customer);
        }
        *it = (denom == 0.0) ? 0.0 : (*it / denom);
    }
    return lambda;
}

arma::mat fastcascade_esri(
    const arma::sp_mat& lambda_d,
    const arma::sp_mat& lambda_u,
    const arma::sp_mat& psi_mat,
    const arma::sp_mat& psup,
    arma::mat h_weights,
    double eps,
    const arma::uvec& essential_cols,
    const arma::uvec& nonessential_cols
) {
    const arma::uword n = lambda_d.n_rows;
    const arma::uword m = psup.n_cols;
    const arma::uword n_def = psi_mat.n_cols;
    const arma::uword n_weights = h_weights.n_cols;

    h_weights.each_row() /= arma::sum(h_weights, 0);
    arma::mat ESRI(n_def, 3 * n_weights, arma::fill::zeros);

    for (arma::uword scenario = 0; scenario < n_def; ++scenario) {
        if ((scenario % 100) == 0) {
            std::cout << "Iteration " << scenario << "\n";
        }

        arma::vec psi_d = arma::ones<arma::vec>(n) - arma::vec(psi_mat.col(scenario));
        arma::vec psi_u = psi_d;
        arma::vec h_t_d = psi_d;
        arma::vec h_t_u = psi_u;
        arma::vec h_tp1_d(n, arma::fill::ones);
        arma::vec h_tp1_u(n, arma::fill::ones);
        arma::mat pi_t(n, m, arma::fill::ones);

        bool crit = true;
        while (crit) {
            arma::sp_mat psup_h = scale_each_row(psup, arma::ones<arma::vec>(n) - h_t_d);
            pi_t = arma::mat(n, m, arma::fill::ones) - arma::mat(arma::trans(lambda_d) * psup_h);

            arma::vec leo_inputs(n, arma::fill::ones);
            if (!essential_cols.is_empty()) {
                leo_inputs = arma::min(pi_t.cols(essential_cols), 1);
            }

            arma::vec lin_inputs(n, arma::fill::ones);
            if (!nonessential_cols.is_empty()) {
                lin_inputs = arma::sum(pi_t.cols(nonessential_cols), 1) -
                             arma::vec(n).fill(nonessential_cols.n_elem) + 1.0;
            }

            h_tp1_d = arma::max(
                arma::zeros<arma::vec>(n),
                arma::min(psi_d, arma::min(lin_inputs, leo_inputs))
            );

            h_tp1_u = arma::max(
                arma::zeros<arma::vec>(n),
                arma::min(psi_u, arma::ones<arma::vec>(n) - lambda_u * (arma::ones<arma::vec>(n) - h_t_u))
            );

            crit = arma::as_scalar(arma::max(arma::max(-h_tp1_d + h_t_d, -h_tp1_u + h_t_u))) > eps;
            h_t_d = h_tp1_d;
            h_t_u = h_tp1_u;
        }

        arma::mat mindownup = arma::ones<arma::mat>(n, 3) -
                              arma::join_horiz(arma::min(h_t_d, h_t_u), h_t_d, h_t_u);
        arma::mat weighted = h_weights.t() * mindownup;
        ESRI.row(scenario) = weighted.as_row();
    }

    return ESRI;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "usage: fastcascade_cli <edges.tsv> <industries.tsv> <essential.tsv> <tol> <metrics.tsv> <scores.tsv>\n";
        return 1;
    }

    const std::string edges_path = argv[1];
    const std::string industries_path = argv[2];
    const std::string essential_path = argv[3];
    const double tol = std::stod(argv[4]);
    const std::string metrics_path = argv[5];
    const std::string scores_path = argv[6];

    arma::uvec p_cons = read_uvec_1based(industries_path);
    arma::uvec essential_flags = read_uvec_raw(essential_path);
    const arma::uword n = p_cons.n_elem;
    const arma::uword m = essential_flags.n_elem;

    arma::sp_mat W = read_edges(edges_path, n);
    arma::vec out_str = arma::vec(arma::mat(arma::sum(W, 1)));
    arma::sp_mat lambda_u = row_normalize(W, out_str);
    arma::sp_mat lambda_d = lambda_d_calc(W, p_cons, essential_flags);
    arma::sp_mat psup = build_psup(p_cons, m);

    arma::umat diag_locations(2, n);
    arma::vec diag_values(n, arma::fill::ones);
    for (arma::uword i = 0; i < n; ++i) {
        diag_locations(0, i) = i;
        diag_locations(1, i) = i;
    }
    arma::sp_mat psi_mat(diag_locations, diag_values, n, n);

    arma::mat h_weights = arma::mat(out_str);
    arma::uvec essential_cols = arma::find(essential_flags == 1);
    arma::uvec nonessential_cols = arma::find(essential_flags == 0);

    const auto start = std::chrono::steady_clock::now();
    arma::mat ESRI = fastcascade_esri(
        lambda_d,
        lambda_u,
        psi_mat,
        psup,
        h_weights,
        tol,
        essential_cols,
        nonessential_cols
    );
    const auto end = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(end - start).count();

    std::ofstream metrics_out(metrics_path);
    if (!metrics_out) {
        throw std::runtime_error("failed to open " + metrics_path);
    }
    metrics_out << "diem_total_s\t" << elapsed << "\n";
    metrics_out << "diem_score_min\t" << ESRI.col(0).min() << "\n";
    metrics_out << "diem_score_max\t" << ESRI.col(0).max() << "\n";
    metrics_out << "diem_downstream_score_min\t" << ESRI.col(1).min() << "\n";
    metrics_out << "diem_downstream_score_max\t" << ESRI.col(1).max() << "\n";
    metrics_out << "diem_upstream_score_min\t" << ESRI.col(2).min() << "\n";
    metrics_out << "diem_upstream_score_max\t" << ESRI.col(2).max() << "\n";

    std::ofstream scores_out(scores_path);
    if (!scores_out) {
        throw std::runtime_error("failed to open " + scores_path);
    }
    for (arma::uword scenario = 0; scenario < ESRI.n_rows; ++scenario) {
        scores_out << ESRI(scenario, 0) << '\t'
                   << ESRI(scenario, 1) << '\t'
                   << ESRI(scenario, 2) << "\n";
    }
    return 0;
}
