#include "allo_extent.hpp"
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <set>
#include "../common/helpers/eigen_helpers.hpp"

namespace ear {

  namespace {
    const int num_vs = 40;
    const double NEG130DBEXP_LIM = 6.5;
    const double NEG130DB_LIM = std::pow(10.0, -NEG130DBEXP_LIM);
    const double pi = boost::math::constants::pi<double>();

    double scaleSize(double v) {
      Eigen::Matrix<double, 5, 1> xp, yp;
      xp << 0.0, 0.2, 0.5, 0.75, 1.0;
      yp << 0.0, 0.3, 1.0, 1.8, 2.8;
      return interp(std::min(v, 1.0), xp, yp);
    }

    bool allEqualColumns(const Eigen::MatrixXd& positions, int c0, int c1) {
      for (Eigen::Index i = 0; i < positions.rows(); i++)
        if (positions(i, c0) != positions(0, c0) ||
            positions(i, c1) != positions(0, c1))
          return false;
      return true;
    }

    bool allEqualColumn(const Eigen::MatrixXd& positions, int c) {
      for (Eigen::Index i = 0; i < positions.rows(); i++)
        if (positions(i, c) != positions(0, c)) return false;
      return true;
    }

    double sEff(const Eigen::MatrixXd& channelPositions, double sx, double sy,
                double sz) {
      if (allEqualColumns(channelPositions, 1, 2)) {
        // Speakers in a left/right line
        return sx;
      } else if (allEqualColumn(channelPositions, 2)) {
        // Speakers in a horizontal plane
        double lo = std::min(sx, sy), hi = std::max(sx, sy);
        return (3.0 / 4.0) * hi + (1.0 / 4.0) * lo;
      } else {
        // Speakers in a cube
        std::array<double, 3> sorted{sx, sy, sz};
        std::sort(sorted.begin(), sorted.end());
        return (6.0 / 9.0) * sorted[2] + (2.0 / 9.0) * sorted[1] +
               (1.0 / 9.0) * sorted[0];
      }
    }

    double pExp(double s_eff) {
      if (s_eff <= 0.5) return 6.0;
      double s_max = 2.8;
      return 6.0 - 4.0 * ((s_eff - 0.5) / (s_max - 0.5));
    }

    double h(double c, double s, double d_bound) {
      (void)c;
      if (d_bound >= 2.0 * s && d_bound >= 0.4) {
        double n = std::max(2.0 * s, 0.4);
        double nnn = n * n * n;
        return std::pow(nnn / (0.16 * 2.0 * s), 1.0 / 3.0);
      } else {
        double a = d_bound / 0.4;
        double b = d_bound / 2.0 * (a * a);
        return std::pow(b, 1.0 / 3.0);
      }
    }

    double dBound(int dim, double xo, double yo, double zo) {
      if (dim == 1) return std::min(xo + 1, 1 - xo);
      if (dim == 2)
        return std::min({xo + 1, 1 - xo, yo + 1, 1 - yo});
      return std::min({xo + 1, 1 - xo, yo + 1, 1 - yo, zo + 1, 1 - zo});
    }

    double mu(int dim, double sx, double sy, double sz, double xo, double yo,
              double zo) {
      double d_bound = dBound(dim, xo, yo, zo);
      if (dim == 1) {
        double n = h(xo, sx, d_bound);
        return n * n * n;
      } else if (dim == 2) {
        double n = h(xo, sx, d_bound) * h(yo, sy, d_bound);
        return std::pow(n, 1.5);
      } else {
        return h(xo, sx, d_bound) * h(yo, sy, d_bound) * h(zo, sz, d_bound);
      }
    }

    int dim(const Eigen::MatrixXd& channelPositions) {
      int d = 0;
      for (int c = 0; c < 3; c++)
        if (!allEqualColumn(channelPositions, c)) d++;
      return d;
    }

    Eigen::VectorXd linspace(double start, double stop, int num) {
      Eigen::VectorXd ret(num);
      double step = (stop - start) / (num - 1);
      for (int i = 0; i < num; i++) ret(i) = start + step * i;
      ret(num - 1) = stop;
      return ret;
    }

    /// weights wx, wy, wz for the virtual source grids
    void calcW(double xo, double yo, double zo, double sx, double sy, double sz,
               const Eigen::VectorXd& xs, const Eigen::VectorXd& ys,
               const Eigen::VectorXd& zs, Eigen::VectorXd& wx,
               Eigen::VectorXd& wy, Eigen::VectorXd& wz) {
      wx = (-((1.5 * (xs.array() - xo) / (2.0 * sx)).pow(4.0)).min(NEG130DBEXP_LIM))
               .unaryExpr([](double t) { return std::pow(10.0, t); });
      wy = (-((1.5 * (ys.array() - yo) / (2.0 * sy)).pow(4.0)).min(NEG130DBEXP_LIM))
               .unaryExpr([](double t) { return std::pow(10.0, t); });
      wz = (-((1.5 * (zs.array() - zo) / sz).pow(4.0)).min(NEG130DBEXP_LIM))
               .unaryExpr([](double t) { return std::pow(10.0, t); })
               .array() *
           (zs.array() * pi * (3.0 / 7.0)).cos();
    }

    /// f = sum over the grid of (g_point * w)^p, with small values zeroed
    Eigen::VectorXd calcF(double p, const Eigen::VectorXd& w,
                          const Eigen::MatrixXd& g_point) {
      Eigen::VectorXd f =
          (g_point.array().rowwise() * w.transpose().array()).pow(p).rowwise().sum();
      return f.unaryExpr([](double v) { return v < NEG130DB_LIM ? 0.0 : v; });
    }

    struct Bounds {
      bool hasLo = false, hasHi = false;
      double lo = 0.0, hi = 0.0;
      void addLo(double v) {
        lo = hasLo ? std::max(lo, v) : v;
        hasLo = true;
      }
      void addHi(double v) {
        hi = hasHi ? std::min(hi, v) : v;
        hasHi = true;
      }
    };

    Bounds findPlaneZ(double z, const Eigen::MatrixXd& positions) {
      Bounds b;
      for (Eigen::Index i = 0; i < positions.rows(); i++) {
        double s_z = positions(i, 2);
        if (s_z <= z) b.addLo(s_z);
        if (s_z >= z) b.addHi(s_z);
      }
      return b;
    }

    Bounds findRowY(double y, double plane_z, const Eigen::MatrixXd& positions) {
      Bounds b;
      for (Eigen::Index i = 0; i < positions.rows(); i++) {
        if (positions(i, 2) == plane_z) {
          double s_y = positions(i, 1);
          if (s_y <= y) b.addLo(s_y);
          if (s_y >= y) b.addHi(s_y);
        }
      }
      return b;
    }

    Bounds findColumnX(double x, double row_y, double plane_z,
                       const Eigen::MatrixXd& positions) {
      Bounds b;
      for (Eigen::Index i = 0; i < positions.rows(); i++) {
        if (positions(i, 2) == plane_z && positions(i, 1) == row_y) {
          double s_x = positions(i, 0);
          if (s_x <= x) b.addLo(s_x);
          if (s_x >= x) b.addHi(s_x);
        }
      }
      return b;
    }

    /// gain for one loudspeaker coordinate given the bounds around the
    /// virtual source coordinate
    double gainForBounds(const Bounds& b, double pos, double v) {
      if (!b.hasLo) return pos != b.hi ? 0.0 : 1.0;
      if (!b.hasHi) return pos != b.lo ? 0.0 : 1.0;
      if (b.lo <= pos && pos <= b.hi) {
        if (b.lo == b.hi) return 1.0;
        if (b.lo == pos)
          return std::cos((v - b.lo) / (b.hi - b.lo) * pi / 2.0);
        return std::sin((v - b.lo) / (b.hi - b.lo) * pi / 2.0);
      }
      return 0.0;
    }

    /// per-axis point source gains for each loudspeaker (rows) at each
    /// virtual source coordinate (columns)
    void calcGPointSeparated(const Eigen::MatrixXd& positions,
                             const Eigen::VectorXd& xs,
                             const Eigen::VectorXd& ys,
                             const Eigen::VectorXd& zs, Eigen::MatrixXd& gx,
                             Eigen::MatrixXd& gy, Eigen::MatrixXd& gz) {
      const Eigen::Index n = positions.rows();
      gx.resize(n, xs.size());
      gy.resize(n, ys.size());
      gz.resize(n, zs.size());

      for (Eigen::Index i = 0; i < n; i++) {
        double pos_x = positions(i, 0), pos_y = positions(i, 1),
               pos_z = positions(i, 2);

        for (Eigen::Index k = 0; k < zs.size(); k++)
          gz(i, k) = gainForBounds(findPlaneZ(zs(k), positions), pos_z, zs(k));

        for (Eigen::Index k = 0; k < ys.size(); k++)
          gy(i, k) = gainForBounds(findRowY(ys(k), pos_z, positions), pos_y, ys(k));

        for (Eigen::Index k = 0; k < xs.size(); k++)
          gx(i, k) = gainForBounds(findColumnX(xs(k), pos_y, pos_z, positions),
                                   pos_x, xs(k));
      }
    }

    Eigen::VectorXd safeNorm(const Eigen::VectorXd& vec) {
      double length = vec.norm();
      if (length > 1e-16) return vec / length;
      return Eigen::VectorXd::Zero(vec.size());
    }
  }  // namespace

  Eigen::VectorXd alloExtentGains(const Eigen::MatrixXd& channelPositions,
                                  const Eigen::Vector3d& position,
                                  double sizeX, double sizeY, double sizeZ) {
    double xo = position(0), yo = position(1), zo = position(2);

    std::set<double> distinctZ;
    for (Eigen::Index i = 0; i < channelPositions.rows(); i++)
      distinctZ.insert(channelPositions(i, 2));
    bool threeLayers = distinctZ.size() >= 3;

    const int Nx = num_vs, Ny = num_vs;
    const int Nz = threeLayers ? num_vs : num_vs / 2;
    Eigen::VectorXd xs = linspace(-1.0, 1.0, Nx);
    Eigen::VectorXd ys = linspace(-1.0, 1.0, Ny);
    Eigen::VectorXd zs;
    if (threeLayers) {
      zs = linspace(-1.0, 1.0, Nz);
    } else {
      zs = linspace(0.0, 1.0, Nz);
      zo = std::max(0.0, zo);
    }

    double sx = std::max(scaleSize(sizeX), 2.0 / (Nx - 1));
    double sy = std::max(scaleSize(sizeY), 2.0 / (Ny - 1));
    double sz = std::max(scaleSize(sizeZ), 2.0 / (Nz - 1));
    double s_eff = sEff(channelPositions, sx, sy, sz);
    double p = pExp(s_eff);
    int d = dim(channelPositions);

    double mu_ = mu(d, sx, sy, sz, xo, yo, zo);
    Eigen::VectorXd wx, wy, wz;
    calcW(xo, yo, zo, sx, sy, sz, xs, ys, zs, wx, wy, wz);

    Eigen::MatrixXd g_point_x, g_point_y, g_point_z;
    calcGPointSeparated(channelPositions, xs, ys, zs, g_point_x, g_point_y,
                        g_point_z);
    Eigen::VectorXd fx = calcF(p, wx, g_point_x);
    Eigen::VectorXd fy = calcF(p, wy, g_point_y);
    Eigen::VectorXd fz = calcF(p, wz, g_point_z);

    Eigen::VectorXd g_inside = fx.array() * fy.array() * fz.array();
    Eigen::VectorXd g_inside_norm = safeNorm(g_inside);

    const Eigen::Index last_x = xs.size() - 1, last_y = ys.size() - 1,
                       last_z = zs.size() - 1;
    Eigen::VectorXd b_floor = (g_point_z.col(0) * wz(0)).array().pow(p);
    Eigen::VectorXd b_ceil = (g_point_z.col(last_z) * wz(last_z)).array().pow(p);
    Eigen::VectorXd b_left = (g_point_x.col(0) * wx(0)).array().pow(p);
    Eigen::VectorXd b_right = (g_point_x.col(last_x) * wx(last_x)).array().pow(p);
    Eigen::VectorXd b_front = (g_point_y.col(0) * wy(0)).array().pow(p);
    Eigen::VectorXd b_back = (g_point_y.col(last_y) * wy(last_y)).array().pow(p);

    Eigen::VectorXd g_bound =
        (b_left.array() * fy.array() * fz.array() +
         b_right.array() * fy.array() * fz.array() +
         fx.array() * b_front.array() * fz.array() +
         fx.array() * b_back.array() * fz.array() +
         fx.array() * fy.array() * b_ceil.array() +
         fx.array() * fy.array() * b_floor.array());

    Eigen::VectorXd g_size =
        (g_bound.array() + (mu_ * g_inside_norm.array())).pow(1.0 / p);
    Eigen::VectorXd g_size_norm = safeNorm(g_size);

    double s_fade = 0.2;
    double alpha, beta;
    if (s_eff < s_fade) {
      alpha = std::cos((s_eff * pi) / (s_fade * 2.0));
      beta = std::sin((s_eff * pi) / (s_fade * 2.0));
    } else {
      alpha = 0.0;
      beta = 1.0;
    }

    Eigen::MatrixXd gpx, gpy, gpz;
    Eigen::VectorXd xo_v(1), yo_v(1), zo_v(1);
    xo_v << xo;
    yo_v << yo;
    zo_v << zo;
    calcGPointSeparated(channelPositions, xo_v, yo_v, zo_v, gpx, gpy, gpz);
    Eigen::VectorXd g_point = gpx.col(0).array() * gpy.col(0).array() * gpz.col(0).array();

    Eigen::VectorXd g_total = alpha * g_point + beta * g_size_norm;
    return safeNorm(g_total);
  }

}  // namespace ear
