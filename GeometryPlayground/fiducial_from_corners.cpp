// Pure C++ (no ROOT/STAR dependency -- compiles standalone with g++),
// driven by real detector-corner coordinates read from fcs_corners.txt
// (produced by get_fcs_geometry.C directly from StFcsDb -- i.e. from the
// STAR conditions DB for whatever calibration era that macro was pointed
// at).
//
// Pipeline:
//   1. Read the 4 corners each of ECAL north/south (det 0,1) and HCAL
//      north/south (det 2,3), in STAR global cm.
//   2. Project HCAL corners onto the ECAL plane (z = min ECAL corner z)
//      via simple radial scaling (x*z_plane/z, y*z_plane/z). ECAL's own
//      corners get the same projection: its front face is tilted, not
//      perpendicular to z (confirmed from real corner data -- the inner
//      column, nearest the beam pipe, sits at larger z than the outer
//      column), so only the corners already at z_ecal (the outer column,
//      by construction -- see z_ecal's definition below) can be used
//      as-is; every other corner must be projected too, same as HCAL's.
//   3. Per side (north/south), report HCAL's own projected footprint as
//      the fiducial rectangle (see the "Legacy-matching method" comment
//      below for why), plus a diagnostic cross-check against the true
//      ECAL/HCAL overlap polygon (Sutherland-Hodgman convex intersection).
//
// NOTE: no buffer/erosion is applied anywhere in this tool -- the printed
// x_inner/x_outer/y_min/y_max are the raw, un-buffered rectangle, matching
// JetParameters.h's cut_x_inner/cut_x_outer/cut_y_min/cut_y_max exactly.
// Production applies its buffer live at cut-check time instead (see
// reco_fiducial_buffer/cut_buffer and pass_fiducial_cut() in
// JetParameters.h) -- do NOT erode here too, or jets would get buffered
// twice.
//
// Usage: ./fiducial_from_corners [corners.txt]

#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <limits>

struct Pt2 { double x, y; };

static double cross2(const Pt2 &O, const Pt2 &A, const Pt2 &B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

static double polygonSignedArea(const std::vector<Pt2> &poly) {
    double a = 0.0;
    size_t n = poly.size();
    for (size_t i = 0; i < n; i++) {
        const Pt2 &p1 = poly[i];
        const Pt2 &p2 = poly[(i + 1) % n];
        a += p1.x * p2.y - p2.x * p1.y;
    }
    return a / 2.0;
}

static std::vector<Pt2> ensureCCW(std::vector<Pt2> poly) {
    if (polygonSignedArea(poly) < 0.0) std::reverse(poly.begin(), poly.end());
    return poly;
}

// Keep the region to the left of directed line A->B (i.e. cross2(A,B,P) >= 0).
static std::vector<Pt2> clipHalfPlane(const std::vector<Pt2> &subject, Pt2 A, Pt2 B) {
    std::vector<Pt2> out;
    size_t n = subject.size();
    if (n == 0) return out;
    for (size_t i = 0; i < n; i++) {
        Pt2 cur = subject[i];
        Pt2 prev = subject[(i + n - 1) % n];
        double curSide = cross2(A, B, cur);
        double prevSide = cross2(A, B, prev);
        bool curIn = curSide >= 0.0;
        bool prevIn = prevSide >= 0.0;
        if (curIn) {
            if (!prevIn) {
                double t = prevSide / (prevSide - curSide);
                out.push_back(Pt2{prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y)});
            }
            out.push_back(cur);
        } else if (prevIn) {
            double t = prevSide / (prevSide - curSide);
            out.push_back(Pt2{prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y)});
        }
    }
    return out;
}

static std::vector<Pt2> intersectConvex(std::vector<Pt2> A, std::vector<Pt2> B) {
    A = ensureCCW(A);
    B = ensureCCW(B);
    std::vector<Pt2> result = A;
    size_t n = B.size();
    for (size_t i = 0; i < n && !result.empty(); i++) {
        result = clipHalfPlane(result, B[i], B[(i + 1) % n]);
    }
    return result;
}

static void boundingBox(const std::vector<Pt2> &poly, double &xmin, double &xmax, double &ymin, double &ymax) {
    xmin = ymin = std::numeric_limits<double>::infinity();
    xmax = ymax = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < poly.size(); i++) {
        xmin = std::min(xmin, poly[i].x);
        xmax = std::max(xmax, poly[i].x);
        ymin = std::min(ymin, poly[i].y);
        ymax = std::max(ymax, poly[i].y);
    }
}

static void printPoly(const char *label, const std::vector<Pt2> &poly) {
    printf("%s (%zu pts):", label, poly.size());
    for (size_t i = 0; i < poly.size(); i++) printf(" (%.4f, %.4f)", poly[i].x, poly[i].y);
    printf("\n");
}

int main(int argc, char **argv) {
    std::string infile = (argc > 1) ? argv[1] : "fcs_corners.txt";

    // corners[det] = 4 (x,y,z) triples, in the fixed col1row1 -> colMax,row1
    // -> colMax,rowMax -> col1,rowMax order get_fcs_geometry.C wrote them in.
    std::vector<Pt2> corners2d[4];
    std::vector<double> cornersZ[4];

    std::ifstream fin(infile.c_str());
    if (!fin) {
        fprintf(stderr, "Error: cannot open %s\n", infile.c_str());
        return 1;
    }
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int det; double x, y, z;
        if (!(iss >> det >> x >> y >> z)) continue;
        if (det < 0 || det > 3) continue;
        corners2d[det].push_back(Pt2{x, y});
        cornersZ[det].push_back(z);
    }
    for (int d = 0; d < 4; d++) {
        if (corners2d[d].size() != 4) {
            fprintf(stderr, "Error: detector %d has %zu corners, expected 4\n", d, corners2d[d].size());
            return 1;
        }
    }

    // ECAL front plane: min z among the two ECAL detectors' corners.
    double z_ecal = std::numeric_limits<double>::infinity();
    for (int d = 0; d <= 1; d++)
        for (size_t i = 0; i < cornersZ[d].size(); i++)
            z_ecal = std::min(z_ecal, cornersZ[d][i]);
    printf("Projection plane z_ecal = %.4f cm\n", z_ecal);

    // Project HCAL (det 2 -> north/left, det 3 -> south/right) onto z_ecal.
    std::vector<Pt2> hcalProj[2]; // [0]=det2, [1]=det3
    for (int side = 0; side < 2; side++) {
        int det = side + 2;
        for (size_t i = 0; i < corners2d[det].size(); i++) {
            double x = corners2d[det][i].x, y = corners2d[det][i].y, z = cornersZ[det][i];
            hcalProj[side].push_back(Pt2{x * z_ecal / z, y * z_ecal / z});
        }
    }

    // Project ECAL's OWN corners onto z_ecal too. ECAL's front face is
    // tilted, not perpendicular to z -- confirmed from real corner data,
    // where the inner-column (nearest beam pipe) corners sit at larger z
    // than the outer-column ones. z_ecal is defined as the closest-to-origin
    // edge of that tilted face (min corner z), so a corner's raw (x,y) is
    // only where the origin-ray through it crosses z_ecal if that corner
    // already IS at z_ecal (true for the outer column here, not the inner
    // one). Every other corner must be radially projected from the origin,
    // same as HCAL's corners above -- this matches how jet positions
    // themselves get projected (z_proj in JetParameters.h).
    std::vector<Pt2> ecalProj[2]; // [0]=det0, [1]=det1
    for (int side = 0; side < 2; side++) {
        int det = side;
        for (size_t i = 0; i < corners2d[det].size(); i++) {
            double x = corners2d[det][i].x, y = corners2d[det][i].y, z = cornersZ[det][i];
            ecalProj[side].push_back(Pt2{x * z_ecal / z, y * z_ecal / z});
        }
    }

    const char *sideNames[2] = {"north (left)", "south (right)"};

    // Legacy-matching method (confirmed against production JetParameters.h
    // to 5+ decimal places at the nominal era): the fiducial rectangle is
    // HCAL's OWN projected footprint, not the true ECAL/HCAL overlap
    // polygon -- ECAL is normally the larger, non-limiting detector, and
    // using HCAL's inner column (nearest the beam pipe, col=1) for the y
    // bound is the conservative/narrower choice, since |y| grows toward
    // the outer column as the projection scale factor increases. x_inner/
    // x_outer are just HCAL's col1/colMax projected x (constant along a
    // column, so row-independent).
    //
    // IMPORTANT: north and south are NOT assumed symmetric -- each side's
    // rectangle is reported and used independently. Combining them into a
    // single reflected rectangle (as this tool originally did) silently
    // discards any north/south asymmetry in the real detector position.
    double side_x_inner[2], side_x_outer[2], side_y_min[2], side_y_max[2];

    // ECAL's own footprint (from ecalProj, corners projected onto z_ecal --
    // see the projection comment above), used for the ECAL-only fiducial
    // boundary (EM-only reco jets, no HCAL). Same col1/colMax-column,
    // min/max-row bounding-box convention as the HCAL-footprint rectangle
    // above.
    double side_x_inner_ecal[2], side_x_outer_ecal[2], side_y_min_ecal[2], side_y_max_ecal[2];

    // corner index order is col1,row1=0 ; colMax,row1=1 ; colMax,rowMax=2 ; col1,rowMax=3
    const int COL1_ROW1 = 0, COLMAX_ROW1 = 1, COLMAX_ROWMAX = 2, COL1_ROWMAX = 3;

    for (int side = 0; side < 2; side++) {
        int ecalDet = side; // 0 or 1
        printf("\n--- %s (ECAL det %d / HCAL det %d) ---\n", sideNames[side], ecalDet, side + 2);
        printPoly("ECAL", corners2d[ecalDet]);
        printPoly("HCAL projected", hcalProj[side]);

        double hx_in  = hcalProj[side][COL1_ROW1].x;
        double hx_out = hcalProj[side][COLMAX_ROW1].x;
        double hy_at_col1_row1  = hcalProj[side][COL1_ROW1].y;
        double hy_at_col1_rowMax = hcalProj[side][COL1_ROWMAX].y;
        double ax_inner = std::min(std::fabs(hx_in), std::fabs(hx_out));
        double ax_outer = std::max(std::fabs(hx_in), std::fabs(hx_out));
        double ymin = std::min(hy_at_col1_row1, hy_at_col1_rowMax);
        double ymax = std::max(hy_at_col1_row1, hy_at_col1_rowMax);
        printf("  HCAL-footprint rectangle: |x| in [%.8f, %.8f], y in [%.8f, %.8f]\n", ax_inner, ax_outer, ymin, ymax);

        // Diagnostic cross-check against the true ECAL/HCAL overlap-polygon
        // bbox. The legacy method (above) relies on an implicit assumption:
        // ECAL's raw rectangle is always wider than HCAL's projected
        // footprint, so HCAL alone determines the fiducial boundary. y
        // is EXPECTED to differ a few mm from the true overlap bbox by
        // design (see header comment) -- that's not a red flag. x is not
        // expected to differ at all: if it does, ECAL has become the
        // limiting detector in x for this era, and the legacy method's
        // underlying assumption no longer holds.
        std::vector<Pt2> overlap = intersectConvex(corners2d[ecalDet], hcalProj[side]);
        double oxmin, oxmax, oymin, oymax;
        boundingBox(overlap, oxmin, oxmax, oymin, oymax);
        double true_ax_inner = std::min(std::fabs(oxmin), std::fabs(oxmax));
        double true_ax_outer = std::max(std::fabs(oxmin), std::fabs(oxmax));
        printf("  [diagnostic] true ECAL/HCAL overlap bbox: |x| in [%.8f, %.8f], y in [%.8f, %.8f]\n",
               true_ax_inner, true_ax_outer, oymin, oymax);
        const double tol = 1e-3; // cm
        if (std::fabs(true_ax_inner - ax_inner) > tol || std::fabs(true_ax_outer - ax_outer) > tol) {
            printf("  *** WARNING: x_inner/x_outer diverge from the true overlap bbox by >%.g cm here.\n"
                   "      ECAL has become the limiting detector in x for this era -- the legacy\n"
                   "      method's assumption (HCAL alone determines the boundary) no longer holds;\n"
                   "      review before trusting the numbers below. ***\n", tol);
        }

        side_x_inner[side] = ax_inner;
        side_x_outer[side] = ax_outer;
        side_y_min[side] = ymin;
        side_y_max[side] = ymax;

        // ECAL's own footprint, from its corners projected onto z_ecal
        // (ecalProj, not the raw corners2d -- see the projection comment
        // above; the outer column happens to already sit at z_ecal so it's
        // unaffected, but the inner column does not and must be projected).
        printPoly("ECAL projected", ecalProj[side]);
        double ex_in  = ecalProj[side][COL1_ROW1].x;
        double ex_out = ecalProj[side][COLMAX_ROW1].x;
        double ey_at_col1_row1   = ecalProj[side][COL1_ROW1].y;
        double ey_at_col1_rowMax = ecalProj[side][COL1_ROWMAX].y;
        double eax_inner = std::min(std::fabs(ex_in), std::fabs(ex_out));
        double eax_outer = std::max(std::fabs(ex_in), std::fabs(ex_out));
        double eymin = std::min(ey_at_col1_row1, ey_at_col1_rowMax);
        double eymax = std::max(ey_at_col1_row1, ey_at_col1_rowMax);
        printf("  ECAL-only rectangle:      |x| in [%.8f, %.8f], y in [%.8f, %.8f]\n", eax_inner, eax_outer, eymin, eymax);

        side_x_inner_ecal[side] = eax_inner;
        side_x_outer_ecal[side] = eax_outer;
        side_y_min_ecal[side] = eymin;
        side_y_max_ecal[side] = eymax;
    }

    printf("\n=== Asymmetric fiducial rectangles, per side (JetParameters.h units, cm) ===\n");
    printf("static const float cut_x_inner_north = %.8ff;\n", side_x_inner[0]);
    printf("static const float cut_x_outer_north = %.8ff;\n", side_x_outer[0]);
    printf("static const float cut_y_min_north   = %.8ff;\n", side_y_min[0]);
    printf("static const float cut_y_max_north   = %.8ff;\n", side_y_max[0]);
    printf("static const float cut_x_inner_south = %.8ff;\n", side_x_inner[1]);
    printf("static const float cut_x_outer_south = %.8ff;\n", side_x_outer[1]);
    printf("static const float cut_y_min_south   = %.8ff;\n", side_y_min[1]);
    printf("static const float cut_y_max_south   = %.8ff;\n", side_y_max[1]);
    printf("(no buffer applied here -- production applies reco_fiducial_buffer/cut_buffer\n");
    printf(" live at cut-check time instead, see pass_fiducial_cut() in JetParameters.h.)\n");

    if (std::fabs(side_x_inner[0] - side_x_inner[1]) > 1e-3 ||
        std::fabs(side_x_outer[0] - side_x_outer[1]) > 1e-3 ||
        std::fabs(side_y_min[0] - side_y_min[1]) > 1e-3 ||
        std::fabs(side_y_max[0] - side_y_max[1]) > 1e-3) {
        printf("\nNote: north and south differ (as expected -- do not average/reflect them).\n");
    }

    printf("\n=== ECAL-only fiducial rectangles, per side (JetParameters.h units, cm) ===\n");
    printf("static const float cut_x_inner_ecal_north = %.8ff;\n", side_x_inner_ecal[0]);
    printf("static const float cut_x_outer_ecal_north = %.8ff;\n", side_x_outer_ecal[0]);
    printf("static const float cut_y_min_ecal_north   = %.8ff;\n", side_y_min_ecal[0]);
    printf("static const float cut_y_max_ecal_north   = %.8ff;\n", side_y_max_ecal[0]);
    printf("static const float cut_x_inner_ecal_south = %.8ff;\n", side_x_inner_ecal[1]);
    printf("static const float cut_x_outer_ecal_south = %.8ff;\n", side_x_outer_ecal[1]);
    printf("static const float cut_y_min_ecal_south   = %.8ff;\n", side_y_min_ecal[1]);
    printf("static const float cut_y_max_ecal_south   = %.8ff;\n", side_y_max_ecal[1]);

    if (std::fabs(side_x_inner_ecal[0] - side_x_inner_ecal[1]) > 1e-3 ||
        std::fabs(side_x_outer_ecal[0] - side_x_outer_ecal[1]) > 1e-3 ||
        std::fabs(side_y_min_ecal[0] - side_y_min_ecal[1]) > 1e-3 ||
        std::fabs(side_y_max_ecal[0] - side_y_max_ecal[1]) > 1e-3) {
        printf("\nNote: ECAL-only north and south differ (as expected -- do not average/reflect them).\n");
    }

    return 0;
}
