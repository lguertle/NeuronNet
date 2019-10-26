#include "network.h"
#include "random.h"
#include <vector>
#include <math.h>

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::pair<size_t, double> Network::degree(const size_t& n) const
{
	std::pair<size_t, double> number_intensity;
	double intensity(0);
	for(auto x:neighbors(n))
	{
		intensity += x.second;
	}
	number_intensity.first = neighbors(n).size();
	number_intensity.second = intensity;

	return number_intensity;
}

std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& n) const
{
	std::vector<std::pair<size_t, double>> n_neighbors;
	for (auto x:links)
	{
		if (x.first.first == n)
		{
			n_neighbors.push_back(std::pair<size_t, double>(x.first.second, x.second));
		}
	}
	return n_neighbors;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

std::set<size_t> Network::step(const std::vector<double>& thalamic_input)
{
	std::set<size_t> firing_n;
	/*
	double size(thalamic_input.size());
	double sum(0);
	double mean(0);
	double var(0);
	double variance(0);
	for(size_t i(0); i<size; ++i)
	{
		sum += thalamic_input[i];
	}
	mean = sum/size;
	
	for(size_t i(0); i<size; ++i)
	{
		var += pow(thalamic_input[i]-mean,2);
	}
	variance = var/size;*/

	double w(0);
	double activ(0);
	double inhib(0);
	for(size_t i(0); i<neurons.size(); ++i)
	{
		neurons[i].step();
		if (neurons[i].is_inhibitory()) {
			w = 2;
		}
		else {
			w = 5;
		}
			for(auto nb:neighbors(i))
			{
				if(neurons[nb.first].firing())
				{
					if (neurons[nb.first].is_inhibitory()) {
						inhib += nb.second;
					}
					else {
						activ += nb.second;
					}
					firing_n.insert(nb.first);
					neurons[nb.first].reset();
				}
			}
		firing_n.insert(i);
		neurons[i].input(thalamic_input[i]/5*w + 0.5*activ-inhib);
	}
	return firing_n;
}


void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}
