#include <cmath>
#include <random>
#include <vector>

int mod(int x, unsigned m) {
	int r = x % static_cast<int>(m);
	if (r < 0) {
		return r + m;
	} else {
		return r;
	}
}

struct Index3 {
		int x;
		int y;
		int z;
};

struct Lattice {
	std::vector<bool> _occupied;
	std::vector<Index3> _prev_site;
	std::vector<Index3> _next_site;
	unsigned N;
	unsigned k;

	Lattice(unsigned N, unsigned k) :
		_occupied(k * N * N * N, false),
		_prev_site(k * N * N * N, Index3 { }),
		_next_site(k * N * N * N, Index3 { }),
	       N(N), k(k) { }

	bool occupied(int t, Index3 index) const {
		t = mod(t, k);
		index.x = mod(index.x, N);
		index.y = mod(index.y, N);
		index.z = mod(index.z, N);
		return _occupied[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	std::vector<bool>::reference occupied(int t, Index3 index) {
		t = mod(t, k);
		index.x = mod(index.x, N);
		index.y = mod(index.y, N);
		index.z = mod(index.z, N);
		return _occupied[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3 const& prev_site(int t, Index3 index) const {
		t = mod(t, k);
		index.x = mod(index.x, N);
		index.y = mod(index.y, N);
		index.z = mod(index.z, N);
		return _prev_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3& prev_site(int t, Index3 index) {
		t = mod(t, k);
		index.x = mod(index.x, N);
		index.y = mod(index.y, N);
		index.z = mod(index.z, N);
		return _prev_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3 const& next_site(int t, Index3 index) const {
		t = mod(t, k);
		index.x = mod(index.x, N);
		index.y = mod(index.y, N);
		index.z = mod(index.z, N);
		return _next_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3& next_site(int t, Index3 index) {
		t = mod(t, k);
		index.x = mod(index.x, N);
		index.y = mod(index.y, N);
		index.z = mod(index.z, N);
		return _next_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
};

int main(int argc, char** argv) {
	// Initialize lattice.
	unsigned N = 2;    // Number of lattice sites along one axis.
	unsigned k = 1090; // Number of time steps.
	unsigned d = 3;    // Dimension of lattice.
	double t = 1.;     // Hopping constant.
	double mu = 1.;    // Chemical potential.
	double beta = 1.;  // Beta.
	double epsilon = beta / k;
	unsigned event_count = 100000;
	Lattice lattice(N, k);
	// Random number engine.
	std::default_random_engine rng;
	std::uniform_real_distribution<double> prob_dist(0., 1.);
	std::uniform_int_distribution<int> N_dist(0, N - 1);
	std::uniform_int_distribution<int> k_dist(0, k - 1);
	std::uniform_int_distribution<unsigned> dir_dist(0, d - 1);
	bool forwards = true;
	Index3 next_xs = { 0, 0, 0 };
	for (unsigned idx = 0; idx < event_count; ++idx) {
		// Pick point.
		Index3 xs = { N_dist(rng), N_dist(rng), N_dist(rng) };
		int t = k_dist(rng);
		// Check if occupied to set our initial direction correctly.
		forwards = !lattice.occupied(t, xs);
		if (lattice.occupied(t, xs)) {
				forwards = false;
				next_xs = lattice.prev_site(t, xs);
		} else {
				forwards = true;
		}
		while (true) {
			// Evolve.
			if (forwards) {
				if (prob_dist(rng) < std::exp(mu * epsilon)) {
					// Accept.
					// Choose direction.
					next_xs = xs;
					if (prob_dist(rng) > 1. - 2. * t * d * epsilon) {
							// Hop to different site.
							unsigned dir = dir_dist(rng);
							int sign = prob_dist(rng) < 0.5 ? 1 : -1;
							switch (dir) {
							case 0:
								next_xs.x += sign;
								break;
							case 1:
								next_xs.y += sign;
								break;
							case 2:
								next_xs.z += sign;
								break;
							}
					}
					bool old_occupied = lattice.occupied(t + 1, next_xs);
					Index3 old_prev_site = lattice.prev_site(t + 1, next_xs);
					// Create a path to next site.
					lattice.occupied(t + 1, next_xs) = true;
					lattice.prev_site(t + 1, next_xs) = xs;
					lattice.next_site(t, xs) = next_xs;
					xs = next_xs;
					t += 1;
					if (old_occupied) {
							// Flip direction and start moving backwards.
							forwards = false;
							next_xs = old_prev_site;
					} else {
							// Good.
					}
				} else {
					// Reject.
					forwards = false;
					next_xs = lattice.prev_site(t, xs);
				}
			} else {
				// Moving backwards.
				if (prob_dist(rng) < std::exp(-mu * epsilon)) {
					// Accept.
					
				} else {
					// Reject.
				}
			}
			// Check finished.
		}
		// Compute observables.
		// Statistics.
	}
	// END.
}

