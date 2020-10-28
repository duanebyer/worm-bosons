#include <cmath>
#include <ctime>
#include <iostream>
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

	bool operator==(Index3 other) const {
		return x == other.x && y == other.y && z == other.z;
	}
	bool operator!=(Index3 other) const {
		return !(*this == other);
	}
	int operator[](unsigned coord) const {
		switch (coord) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		}
	}
	int& operator[](unsigned coord) {
		switch (coord) {
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		}
	}
};

Index3 mod(Index3 xs, unsigned m) {
	xs.x = mod(xs.x, m);
	xs.y = mod(xs.y, m);
	xs.y = mod(xs.y, m);
	return xs;
}

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

	std::vector<bool>::const_reference occupied(int t, Index3 index) const {
		t = mod(t, k);
		index = mod(index, N);
		return _occupied[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	std::vector<bool>::reference occupied(int t, Index3 index) {
		t = mod(t, k);
		index = mod(index, N);
		return _occupied[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3 const& prev_site(int t, Index3 index) const {
		t = mod(t, k);
		index = mod(index, N);
		return _prev_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3& prev_site(int t, Index3 index) {
		t = mod(t, k);
		index = mod(index, N);
		return _prev_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3 const& next_site(int t, Index3 index) const {
		t = mod(t, k);
		index = mod(index, N);
		return _next_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Index3& next_site(int t, Index3 index) {
		t = mod(t, k);
		index = mod(index, N);
		return _next_site[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
};

int main(int argc, char** argv) {
	// Initialize lattice.
	unsigned N = 2;    // Number of lattice sites along one axis.
	unsigned K = 1090; // Number of time steps.
	unsigned D = 3;    // Dimension of lattice.
	double T = 1.;     // Hopping constant.
	double mu = 1.4;    // Chemical potential.
	double beta = 12.;  // Beta.
	double epsilon = beta / K;
	unsigned event_count = 1000;
	Lattice lattice(N, K);
	// Random number engine.
	std::default_random_engine rng;
	rng.seed(std::time(0));
	std::uniform_real_distribution<double> prob_dist(0., 1.);
	std::uniform_int_distribution<int> N_dist(0, N - 1);
	std::uniform_int_distribution<int> K_dist(0, K - 1);
	std::uniform_int_distribution<unsigned> dir_dist(0, D - 1);
	// Statistics.
	std::vector<double> energies;
	std::vector<unsigned> numbers;
	std::vector<unsigned> winding_numbers;
	// State variables.
	bool forwards = true;
	Index3 next_xs = { 0, 0, 0 };
	Index3 winding = { 0, 0, 0 };
	unsigned hop_space = 0;
	for (unsigned idx = 0; idx < event_count; ++idx) {
		// Pick point.
		Index3 start_xs = { N_dist(rng), N_dist(rng), N_dist(rng) };
		int start_t = K_dist(rng);
		Index3 xs = start_xs;
		int t = start_t;
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
					if (prob_dist(rng) > 1. - 2. * T * D * epsilon) {
						// Hop to different site.
						hop_space += 1;
						unsigned dir = dir_dist(rng);
						int sign = prob_dist(rng) < 0.5 ? 1 : -1;
						next_xs[dir] += sign;
						if (next_xs[dir] < 0) {
							next_xs[dir] += N;
							winding[dir] += 1;
						}
						if (next_xs[dir] >= N) {
							next_xs[dir] -= N;
							winding[dir] -= 1;
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
					if (xs != next_xs) {
						hop_space -= 1;
					}
					lattice.occupied(t, xs) = false;
					xs = next_xs;
					t -= 1;
					next_xs = lattice.prev_site(t, xs);
				} else {
					// Reject.
					forwards = true;
					next_xs = xs;
					if (prob_dist(rng) > 1. - 2. * T * D * epsilon) {
						// Hop to different site.
						hop_space += 1;
						unsigned dir = dir_dist(rng);
						int sign = prob_dist(rng) < 0.5 ? 1 : -1;
						next_xs[dir] += sign;
						if (next_xs[dir] < 0) {
							next_xs[dir] += N;
							winding[dir] += 1;
						}
						if (next_xs[dir] >= N) {
							next_xs[dir] -= N;
							winding[dir] -= 1;
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
				}
			}
			// Check finished.
			if (start_t == mod(t, K) && start_xs == mod(xs, N)) {
				break;
			}
		}
		// Compute observables.
		// Particle number.
		unsigned n = 0;
		for (int z = 0; z < N; ++z) {
			for (int y = 0; y < N; ++y) {
				for (int x = 0; x < N; ++x) {
					if (lattice.occupied(0, Index3 { x, y, z })) {
						n += 1;
					}
				}
			}
		}
		double energy = -hop_space / beta + 2. * T * D * n;
		// Statistics.
		energies.push_back(energy);
		numbers.push_back(n);
		winding_numbers.push_back(winding.x + winding.y + winding.z);
	}
	// END.
	// Compute means.
	double mean_energy = 0.;
	for (double energy : energies) {
		mean_energy += energy;
	}
	mean_energy /= energies.size();
	double mean_number = 0;
	for (unsigned number : numbers) {
		mean_number += number;
	}
	mean_number /= numbers.size();

	std::cout << "Mean energy: " << mean_energy << std::endl;
	std::cout << "Mean number: " << mean_number << std::endl;

	return 0;
}

