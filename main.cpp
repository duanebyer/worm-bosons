#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <string>
#include <vector>

bool const LOG = false;

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
							}
						}
						switch (occupied(t, { x, y, z })) {
						case Occupied::EMPTY:
							std::cout << '.';
							break;
						case Occupied::LINE:
						case Occupied::TAIL:
						case Occupied::MERGE:
							if (mod(next.x, N) == x + 1) {
								std::cout << "→";
							} else if(mod(next.x, N) == x - 1) {
								std::cout << "←";
							} else if (mod(next.y, N) == y + 1) {
								std::cout << "↓";
							} else if (mod(next.y, N) == y - 1) {
								std::cout << "↑";
							} else if (mod(next.z, N) == z + 1) {
								std::cout << "⊙";
							} else if (mod(next.z, N) == z - 1) {
								std::cout << "⊗";
							}
							break;
						case Occupied::HEAD:
							std::cout << 'H';
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
	unsigned K = static_cast<unsigned>(std::ceil(beta / epsilon));
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
	std::vector<unsigned> numbers;
	std::vector<unsigned> winding_numbers;
	// State variables.
	bool forwards = true;
	Index3 merge_prev_xs = { 0, 0, 0 };
	Index3 winding = { 0, 0, 0 };
	int hop_space = 0;
	for (unsigned idx = 0; idx < event_count + burn_count; ++idx) {
		if (idx % (event_count / 100) == 0) {
			std::cout << "Produced " << 100 * static_cast<double>(idx) / event_count << "%" << std::endl;
		}
		if (LOG) {
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
			if (LOG) {
				std::cout << "hops: " << hop_space << std::endl;
				std::cout << "pos: " << t << ", (" << xs[0] << ", " << xs[1] << ", " << xs[2] << ")" << std::endl;
				lattice.print_lattice();
			}
			// Evolve.
			if (forwards) {
				if (prob_dist(rng) <= std::exp(mu * epsilon)) {
					if (LOG) {
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
								winding[dir] += 1;
							}
							if (next_xs[dir] >= N) {
								next_xs[dir] -= N;
								winding[dir] -= 1;
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
							if (LOG) {
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
					if (LOG) {
						std::cout << "forwards reject" << std::endl;
					}
					forwards = false;
					Occupied occupied = lattice.occupied(t, xs);
					if (occupied == Occupied::HEAD) {
					} else if (occupied == Occupied::EMPTY) {
						// Terminate.
						if (LOG) {
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
					if (LOG) {
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
						if (LOG) {
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
					xs = next_xs;
					t -= 1;
					if (t < 0) {
						t += K;
					}
				} else {
					if (LOG) {
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
							winding[dir] += 1;
						}
						if (next_xs[dir] >= N) {
							next_xs[dir] -= N;
							winding[dir] -= 1;
						}
					}
					Occupied next_occupied = lattice.occupied(t, next_xs);
					if (occupied == Occupied::LINE
							|| occupied == Occupied::HEAD
							|| occupied == Occupied::MERGE) {
						if (next_xs == xs) {
							if (occupied == Occupied::LINE) {
								// Terminate.
								if (LOG) {
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
								if (LOG) {
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
		}
		if (LOG) {
			std::cout << "final: " << std::endl;
			lattice.print_lattice();
		}
		// Compute observables.
		// Check validity of lattice.
		if (LOG) {
			for (int t = 0; t < K; ++t) {
				for (int z = 0; z < N; ++z) {
					for (int y = 0; y < N; ++y) {
						for (int x = 0; x < N; ++x) {
							Occupied occupied = lattice.occupied(t, { x, y, z });
							if (occupied != Occupied::LINE && occupied != Occupied::EMPTY) {
								std::cerr
									<< "Error 8: Invalid occupied # "
									<< static_cast<int>(occupied) << " at "
									<< "(" << x << ", " << y << ", " << z << ", " << t << ")" << std::endl;
								return 8;
							}
							if (occupied == Occupied::LINE) {
								Index3 next = lattice.next_site(t, Index3 { x, y, z });
								Occupied next_occupied = lattice.occupied(t + 1, next);
								if (next_occupied == Occupied::EMPTY) {
									std::cerr
										<< "Error 8: Connected from line to empty at "
										<< "(" << x << ", " << y << ", " << z << ", " << t << ")" << std::endl;
									return 8;
								}
								if (mod(lattice.prev_site(t + 1, next), N) != Index3 { x, y, z }) {
									std::cerr
										<< "Error 8: One way connection from "
										<< "(" << x << ", " << y << ", " << z << ", " << t << ")" << std::endl;
									return 8;
								}
							}
						}
					}
				}
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
		// Statistics.
		// Burn the first few configurations.
		if (idx >= burn_count) {
			energies.push_back(energy);
			numbers.push_back(n);
			winding_numbers.push_back(winding.x + winding.y + winding.z);
		}
	}
	// Compute means.
	double mean_energy = 0.;
	for (double energy : energies) {
		mean_energy += energy;
	}
	mean_energy /= event_count;
	double mean_number = 0.;
	for (unsigned number : numbers) {
		mean_number += (double) number;
	}
	mean_number /= event_count;

	std::cout << "Mean energy: " << mean_energy << std::endl;
	std::cout << "Mean number: " << mean_number << std::endl;

	return 0;
}

