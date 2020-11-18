#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

bool const LOG_PROCESS = false;
bool const LOG_STATE = false;
bool const LOG_OBSERVABLES = false;
bool const CHECK_VALIDITY = false;

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
	Index3& operator+=(Index3 const& rhs) {
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	Index3& operator-=(Index3 const& rhs) {
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
};

Index3 operator+(Index3 lhs, Index3 const& rhs) {
	lhs += rhs;
	return lhs;
}

Index3 operator-(Index3 lhs, Index3 const& rhs) {
	lhs -= rhs;
	return lhs;
}

Index3 mod(Index3 xs, unsigned m) {
	xs.x = mod(xs.x, m);
	xs.y = mod(xs.y, m);
	xs.y = mod(xs.y, m);
	return xs;
}

Index3 check_winding(Index3 xs, Index3 next_xs, unsigned N) {
	Index3 winding = { 0, 0, 0 };
	for (unsigned dim = 0; dim < 3; ++dim) {
		if (xs[dim] == 0 && mod(next_xs[dim], N) == N - 1) {
			winding[dim] += 1;
		}
		if (xs[dim] == N - 1 && mod(next_xs[dim], N) == 0) {
			winding[dim] -= 1;
		}
	}
	return winding;
}

// Some helpful functions for computing statistics.
double average(double const* data, std::size_t n) {
	double mean = 0.;
	for (std::size_t idx = 0; idx < n; ++idx) {
		mean += data[idx];
	}
	mean /= n;
	return mean;
}

double variance(double const* data, std::size_t n) {
	double mean = average(data, n);
	double var = 0.;
	for (std::size_t idx = 0; idx < n; ++idx) {
		var += (data[idx] - mean) * (data[idx] - mean);
	}
	var /= n - 1;
	return var;
}

std::vector<double> chunks(double const* data, std::size_t n, std::size_t k) {
	std::size_t pos = 0;
	std::vector<double> chunk;
	while (pos + k <= n) {
		chunk.push_back(average(data + pos, k));
		pos += k;
	}
	return chunk;
}

double correlation(double const* data, std::size_t n, std::size_t k) {
	double corr = 0.;
	double mean = average(data, n);
	for (std::size_t idx = 0; idx < n; ++idx) {
		corr += (data[idx] - mean) * (data[(idx + k) % n] - mean);
	}
	corr /= n;
	return corr;
}

struct EstErr {
	std::size_t chunk_n;
	double est;
	double err;
};

EstErr estimate(std::vector<double> data, double prec = 0.05) {
	double var = variance(data.data(), data.size());
	std::size_t chunk_n = 1;
	// Chunk until correlation has vanished.
	while (std::abs(correlation(data.data(), data.size(), 1)) > prec * var) {
		data = chunks(data.data(), data.size(), 2);
		chunk_n *= 2;
		if (data.size() <= 1) {
			throw std::runtime_error("Not enough data to estimate observables");
		}
	}
	// Return mean and variance of chunks.
	return EstErr {
		chunk_n,
		average(data.data(), data.size()),
		std::sqrt(variance(data.data(), data.size()) / data.size()),
	};
}

enum class Occupied {
	EMPTY,
	LINE,
	HEAD,
	TAIL,
	MERGE,
};

struct Lattice {
	std::vector<Occupied> _occupied;
	std::vector<Index3> _prev_site;
	std::vector<Index3> _next_site;
	unsigned N;
	unsigned k;

	Lattice(unsigned N, unsigned k) :
		_occupied(k * N * N * N, Occupied::EMPTY),
		_prev_site(k * N * N * N, Index3 { }),
		_next_site(k * N * N * N, Index3 { }),
		N(N), k(k) {
	}

	Occupied const& occupied(int t, Index3 index) const {
		t = mod(t, k);
		index = mod(index, N);
		return _occupied[N * N * N * t + N * N * index.x + N * index.y + index.z];
	}
	Occupied& occupied(int t, Index3 index) {
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

	void print_lattice(bool show_occupations=false) const {
		for (unsigned t = 0; t < k; ++t) {
			for (unsigned y = 0; y < N; ++y) {
				for (unsigned z = 0; z < N; ++z) {
					for (unsigned x = 0; x < N; ++x) {
						Index3 next = next_site(t, { x, y, z });
						if (show_occupations) {
							switch (occupied(t, { x, y, z })) {
							case Occupied::EMPTY:
								std::cout << 'E';
								break;
							case Occupied::LINE:
								std::cout << 'L';
								break;
							case Occupied::TAIL:
								std::cout << 'T';
								break;
							case Occupied::HEAD:
								std::cout << 'H';
								break;
							case Occupied::MERGE:
								std::cout << 'M';
								break;
							default:
								std::cout << 'X';
								break;
							}
						}
						switch (occupied(t, { x, y, z })) {
						case Occupied::EMPTY:
							std::cout << '.';
							break;
						case Occupied::LINE:
						case Occupied::TAIL:
						case Occupied::MERGE:
							if (mod(next.x, N) == mod(x + 1, N)) {
								std::cout << "→";
							} else if(mod(next.x, N) == mod(x - 1, N)) {
								std::cout << "←";
							} else if (mod(next.y, N) == mod(y + 1, N)) {
								std::cout << "↓";
							} else if (mod(next.y, N) == mod(y - 1, N)) {
								std::cout << "↑";
							} else if (mod(next.z, N) == mod(z + 1, N)) {
								std::cout << "⊙";
							} else if (mod(next.z, N) == mod(z - 1, N)) {
								std::cout << "⊗";
							} else if (mod(next, N) == Index3 { x, y, z }) {
								std::cout << '*';
							} else {
								std::cout << 'X';
							}
							break;
						case Occupied::HEAD:
							std::cout << 'H';
							break;
						default:
							std::cout << 'X';
							break;
						}
					}
					std::cout << ' ';
				}
				std::cout << std::endl;
			}
			std::cout << "~~~~~~~~" << std::endl;
		}
	}
};

int main(int argc, char** argv) {
	// Command line arguments.
	if (argc != 8) {
		std::cerr << "Usage: worms <# of events> <# to burn> <N> <ε> <μ> <β> <T>" << std::endl;
		return 20;
	}
	unsigned event_count = std::stoi(argv[1]); // # of events.
	unsigned burn_count = std::stoi(argv[2]); // # of events.
	unsigned N = std::stoi(argv[3]);     // Number of lattice sites along one axis.
	unsigned D = 3;                      // Dimension of lattice.
	double epsilon = std::stod(argv[4]); // Time step.
	double mu = std::stod(argv[5]);      // Chemical potential.
	double beta = std::stod(argv[6]);    // Imaginary time.
	double T = std::stod(argv[7]);       // Hopping constant.
	// Number of time steps.
	unsigned K = static_cast<unsigned>(std::round(std::abs(beta) / std::abs(epsilon)));
	if (K <= 2) {
		K = 3;
	}
	epsilon = beta / K;
	// Initialize lattice.
	Lattice lattice(N, K);
	// Random number engine.
	std::default_random_engine rng;
	int seed = std::time(0);
	rng.seed(seed);
	std::cout << "Using random seed " << seed << "." << std::endl;
	std::uniform_real_distribution<double> prob_dist(0., 1.);
	std::uniform_int_distribution<int> N_dist(0, N - 1);
	std::uniform_int_distribution<int> K_dist(0, K - 1);
	std::uniform_int_distribution<unsigned> dir_dist(0, D - 1);
	// Statistics.
	std::vector<double> energies;
	std::vector<double> numbers;
	std::vector<Index3> windings;
	std::vector<double> susceptibilities;
	// State variables.
	bool forwards = true;
	Index3 merge_prev_xs = { 0, 0, 0 };
	Index3 winding = { 0, 0, 0 };
	int hop_space = 0;
	std::cout << "Burning" << std::endl;
	for (unsigned idx = 0; idx < event_count + burn_count; ++idx) {
		if (idx >= burn_count && (idx - burn_count) % (event_count / 100) == 0) {
			std::cout << "Produced " << static_cast<int>(100 * static_cast<double>(idx - burn_count) / event_count) << "%\r";
			std::cout.flush();
		}
		if (LOG_PROCESS || LOG_STATE || LOG_OBSERVABLES) {
			std::cout << "NEW EVENT" << std::endl;
			std::cout << "---------" << std::endl;
		}
		// Pick point.
		Index3 start_xs = { N_dist(rng), N_dist(rng), N_dist(rng) };
		int start_t = K_dist(rng);
		Index3 xs = start_xs;
		int t = start_t;
		// Check if occupied to set our initial direction correctly.
		switch (lattice.occupied(t, xs)) {
		case Occupied::EMPTY:
			forwards = true;
			break;
		case Occupied::LINE:
			forwards = false;
			break;
		default:
			std::cerr << "Error 9" << std::endl;
			return 9;
		}
		bool terminate = false;
		while (!terminate) {
			// Evolve.
			if (forwards) {
				if (prob_dist(rng) <= std::exp(mu * epsilon)) {
					if (LOG_PROCESS) {
						std::cout << "forwards accept" << std::endl;
					}
					// Accept.
					Occupied occupied = lattice.occupied(t, xs);
					if (occupied == Occupied::EMPTY || occupied == Occupied::HEAD) {
						// Choose direction.
						Index3 next_xs = xs;
						if (prob_dist(rng) > 1. - 2. * T * D * epsilon) {
							// Hop to different site.
							hop_space += 1;
							unsigned dir = dir_dist(rng);
							int sign = prob_dist(rng) < 0.5 ? 1 : -1;
							next_xs[dir] += sign;
							if (next_xs[dir] < 0) {
								next_xs[dir] += N;
							}
							if (next_xs[dir] >= N) {
								next_xs[dir] -= N;
							}
						}
						// Create a path to next site.
						Occupied next_occupied = lattice.occupied(t + 1, next_xs);
						if (next_occupied == Occupied::EMPTY) {
							// Propogate.
							lattice.occupied(t + 1, next_xs) = Occupied::HEAD;
						} else if (next_occupied == Occupied::LINE) {
							// Merge.
							forwards = false;
							lattice.occupied(t + 1, next_xs) = Occupied::MERGE;
							merge_prev_xs = lattice.prev_site(t + 1, next_xs);
						} else if (next_occupied == Occupied::TAIL) {
							// Terminate.
							lattice.occupied(t + 1, next_xs) = Occupied::LINE;
							if (LOG_PROCESS) {
								std::cout << "terminate forward hit tail" << std::endl;
							}
							terminate = true;
						} else {
							std::cerr << "Error 1" << std::endl;
							return 1;
						}
						if (occupied == Occupied::EMPTY) {
							lattice.occupied(t, xs) = Occupied::TAIL;
						} else {
							lattice.occupied(t, xs) = Occupied::LINE;
						}
						winding += check_winding(xs, next_xs, N);
						lattice.next_site(t, xs) = next_xs;
						lattice.prev_site(t + 1, next_xs) = xs;
						xs = next_xs;
						t += 1;
						if (t >= K) {
							t -= K;
						}
					} else {
						std::cerr << "Error 2" << std::endl;
						return 2;
					}
				} else {
					if (LOG_PROCESS) {
						std::cout << "forwards reject" << std::endl;
					}
					forwards = false;
					Occupied occupied = lattice.occupied(t, xs);
					if (occupied == Occupied::HEAD) {
					} else if (occupied == Occupied::EMPTY) {
						// Terminate.
						if (LOG_PROCESS) {
							std::cout << "terminate forward immediate" << std::endl;
						}
						terminate = true;
					} else {
						std::cerr << "Error 3" << std::endl;
						return 3;
					}
				}
			} else {
				// Moving backwards.
				if (prob_dist(rng) <= std::exp(-mu * epsilon)) {
					if (LOG_PROCESS) {
						std::cout << "backwards accept" << std::endl;
					}
					// Accept.
					Occupied occupied = lattice.occupied(t, xs);
					Index3 next_xs;
					if (occupied == Occupied::LINE) {
						next_xs = lattice.prev_site(t, xs);
						lattice.occupied(t, xs) = Occupied::TAIL;
					} else if (occupied == Occupied::HEAD) {
						next_xs = lattice.prev_site(t, xs);
						lattice.occupied(t, xs) = Occupied::EMPTY;
					} else if (occupied == Occupied::MERGE) {
						next_xs = merge_prev_xs;
						lattice.occupied(t, xs) = Occupied::LINE;
					} else {
						std::cerr << "Error 4" << std::endl;
						return 4;
					}
					Occupied next_occupied = lattice.occupied(t - 1, next_xs);
					if (next_occupied == Occupied::LINE) {
						lattice.occupied(t - 1, next_xs) = Occupied::HEAD;
					} else if (next_occupied == Occupied::TAIL) {
						lattice.occupied(t - 1, next_xs) = Occupied::EMPTY;
						if (LOG_PROCESS) {
							std::cout << "terminate backwards tail" << std::endl;
						}
						terminate = true;
					} else {
						std::cerr << "Error 5" << std::endl;
						return 5;
					}
					if (xs != next_xs) {
						hop_space -= 1;
					}
					winding -= check_winding(next_xs, xs, N);
					xs = next_xs;
					t -= 1;
					if (t < 0) {
						t += K;
					}
				} else {
					if (LOG_PROCESS) {
						std::cout << "backwards reject" << std::endl;
					}
					// Reject.
					forwards = true;
					Occupied occupied = lattice.occupied(t, xs);
					// Choose direction.
					Index3 next_xs;
					if (occupied == Occupied::LINE || occupied == Occupied::HEAD) {
						next_xs = lattice.prev_site(t, xs);
					} else if (occupied == Occupied::MERGE) {
						next_xs = merge_prev_xs;
					}
					if (xs != next_xs) {
						hop_space -= 1;
					}
					if (prob_dist(rng) > 1. - 2. * T * D * epsilon) {
						// Hop to different site.
						hop_space += 1;
						unsigned dir = dir_dist(rng);
						int sign = prob_dist(rng) < 0.5 ? 1 : -1;
						next_xs[dir] += sign;
						if (next_xs[dir] < 0) {
							next_xs[dir] += N;
						}
						if (next_xs[dir] >= N) {
							next_xs[dir] -= N;
						}
					}
					Occupied next_occupied = lattice.occupied(t, next_xs);
					if (occupied == Occupied::LINE
							|| occupied == Occupied::HEAD
							|| occupied == Occupied::MERGE) {
						if (next_xs == xs) {
							if (occupied == Occupied::LINE) {
								// Terminate.
								if (LOG_PROCESS) {
									std::cout << "terminate backwards immediate" << std::endl;
								}
								terminate = true;
							} else if (occupied == Occupied::HEAD) {
							} else if (occupied == Occupied::MERGE) {
								forwards = false;
								Index3 prev_xs = merge_prev_xs;
								merge_prev_xs = lattice.prev_site(t, next_xs);
								lattice.next_site(t - 1, prev_xs) = next_xs;
								lattice.prev_site(t, next_xs) = prev_xs;
							} else {
								std::cerr << "Error 10" << std::endl;
								return 10;
							}
						} else {
							Index3 prev_xs;
							if (occupied == Occupied::MERGE) {
								prev_xs = merge_prev_xs;
							} else {
								prev_xs = lattice.prev_site(t, xs);
							}
							if (next_occupied == Occupied::EMPTY
									|| next_occupied == Occupied::HEAD) {
								lattice.occupied(t, next_xs) = Occupied::HEAD;
							} else if (next_occupied == Occupied::LINE) {
								forwards = false;
								lattice.occupied(t, next_xs) = Occupied::MERGE;
								merge_prev_xs = lattice.prev_site(t, next_xs);
							} else if ((occupied == Occupied::HEAD || occupied == Occupied::MERGE)
									&& next_occupied == Occupied::TAIL) {
								// Terminate.
								lattice.occupied(t, next_xs) = Occupied::LINE;
								if (LOG_PROCESS) {
									std::cout << "terminate backwards reject" << std::endl;
								}
								terminate = true;
							} else if (next_occupied == Occupied::MERGE) {
								forwards = false;
								merge_prev_xs = lattice.prev_site(t, next_xs);
							} else {
								std::cerr << "Error 6" << std::endl;
								return 6;
							}
							if (occupied == Occupied::LINE) {
								lattice.occupied(t, xs) = Occupied::TAIL;
							} else if (occupied == Occupied::MERGE) {
								lattice.occupied(t, xs) = Occupied::LINE;
							} else {
								lattice.occupied(t, xs) = Occupied::EMPTY;
							}
							winding -= check_winding(prev_xs, xs, N);
							winding += check_winding(prev_xs, next_xs, N);
							lattice.next_site(t - 1, prev_xs) = next_xs;
							lattice.prev_site(t, next_xs) = prev_xs;
							xs = next_xs;
						}
					} else {
						std::cerr << "Error 7" << std::endl;
						return 7;
					}
				}
			}
			if (LOG_PROCESS) {
				std::cout << "hops: " << hop_space << std::endl;
				std::cout << "pos: " << t << ", (" << xs[0] << ", " << xs[1] << ", " << xs[2] << ")" << std::endl;
				std::cout << "wind: " << winding.x << ", " << winding.y << ", " << winding.z << std::endl;
			}
			if (LOG_STATE) {
				lattice.print_lattice();
			}
		}
		// Compute observables.
		// Check validity of lattice.
		unsigned hop_space_validate = 0;
		Index3 winding_validate = { 0, 0, 0 };
		if (CHECK_VALIDITY) {
			for (int t = 0; t < K; ++t) {
				for (int z = 0; z < N; ++z) {
					for (int y = 0; y < N; ++y) {
						for (int x = 0; x < N; ++x) {
							Index3 xs = { x, y, z };
							Occupied occupied = lattice.occupied(t, xs);
							if (occupied != Occupied::LINE && occupied != Occupied::EMPTY) {
								std::cerr
									<< "Error 8: Invalid occupied # "
									<< static_cast<int>(occupied) << " at "
									<< "(" << t << ", " << x << ", " << y << ", " << z << ")" << std::endl;
							}
							if (occupied == Occupied::LINE) {
								Index3 next = lattice.next_site(t, xs);
								Occupied next_occupied = lattice.occupied(t + 1, next);
								if (next_occupied == Occupied::EMPTY) {
									std::cerr
										<< "Error 8: Connected from line to empty at "
										<< "(" << t << ", " << x << ", " << y << ", " << z << ")" << std::endl;
									return 8;
								}
								if (mod(lattice.prev_site(t + 1, next), N) != xs) {
									std::cerr
										<< "Error 8: One way connection from "
										<< "(" << t << ", " << x << ", " << y << ", " << z << ")" << std::endl;
									return 8;
								}
							}
							if (occupied == Occupied::LINE) {
								Index3 next_xs = lattice.next_site(t, xs);
								winding_validate += check_winding(xs, next_xs, N);
								if (mod(next_xs, N) != xs) {
									hop_space_validate += 1;
								}
							}
						}
					}
				}
			}
			if (hop_space != hop_space_validate) {
				std::cerr << "Error 8: Hop space is "
					<< hop_space << "; should be "
					<< hop_space_validate << std::endl;
				return 8;
			}
			if (winding_validate != winding) {
				std::cerr << "Error 8: Winding is "
					<< winding.x << ", " << winding.y << ", " << winding.z
					<< "; should be "
					<< winding_validate.x << ", " << winding_validate.y << ", " << winding_validate.z << std::endl;
				return 8;
			}
		}
		// Particle number.
		unsigned n = 0;
		for (int z = 0; z < N; ++z) {
			for (int y = 0; y < N; ++y) {
				for (int x = 0; x < N; ++x) {
					Occupied occupied = lattice.occupied(0, Index3 { x, y, z });
					if (occupied == Occupied::LINE) {
						n += 1;
					}
				}
			}
		}
		double energy = -hop_space / beta + 2. * T * D * n;
		double susceptibility = 1. / (D * N * beta) * (
			winding.x * winding.x
			+ winding.y * winding.y
			+ winding.z * winding.z);
		// Statistics.
		// Burn the first few configurations.
		if (idx >= burn_count) {
			energies.push_back(energy);
			numbers.push_back(n);
			windings.push_back(winding);
			susceptibilities.push_back(susceptibility);
			if (LOG_OBSERVABLES) {
				std::cout << "e = " << energy << std::endl;
				std::cout << "n = " << n << std::endl;
				std::cout << "ω = " << winding.x << ", " << winding.y << ", " << winding.z << std::endl;
				std::cout << "χ = " << susceptibility << std::endl;
			}
		}
	}
	// Compute means.
	EstErr est_energy = estimate(energies);
	EstErr est_number = estimate(numbers);
	EstErr est_susceptibility = estimate(susceptibilities);

	std::cout << std::endl;
	std::cout << std::scientific << std::setprecision(17);
	std::cout << "Chunk sizes: " << est_energy.chunk_n << ", " << est_number.chunk_n << ", " << est_susceptibility.chunk_n << std::endl;
	std::cout << "Energy: " << est_energy.est << " ± " << est_energy.err << std::endl;
	std::cout << "Number: " << est_number.est << " ± " << est_number.err << std::endl;
	std::cout << "Susceptibility: " << est_susceptibility.est << " ± " << est_susceptibility.err << std::endl;

	return 0;
}

