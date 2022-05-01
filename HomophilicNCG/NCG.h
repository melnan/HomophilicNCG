#pragma once
#include "graph.h"
#include "definitions.h"





template<typename Graph_type, typename cost_f_type>
class Game {
public:
	//constructors
	explicit Game(const Graph_type& input_graph, cost_f_type _alpha, std::set<char>& _list_of_prohib_moves, bool _bilateral) : graph(input_graph), alpha(_alpha), radius_of_steps(graph.graph_num_of_nodes()), list_of_prohib_moves(_list_of_prohib_moves), bilateral(_bilateral)
	{
		infty = static_cast<cost_f_type>(this->graph.infty_value());
	};

	explicit Game(const Graph_type& input_graph, cost_f_type _alpha, int _radius_of_steps, std::set<char>& _list_of_prohib_moves, bool _bilateral) : graph(input_graph), alpha(_alpha), radius_of_steps(_radius_of_steps), list_of_prohib_moves(_list_of_prohib_moves), bilateral(_bilateral)
	{
		infty = static_cast<cost_f_type>(this->graph.infty_value());
	};

	explicit Game(const Graph_type& input_graph, cost_f_type _alpha, bool _bilateral) : graph(input_graph), alpha(_alpha), bilateral(_bilateral) {
		infty = static_cast<cost_f_type>(this->graph.infty_value());
		radius_of_steps = this->graph.graph_num_of_nodes();
		list_of_prohib_moves;
	};

	explicit Game(const Graph_type& input_graph, cost_f_type _alpha, int _radius_of_steps, bool _bilateral) : graph(input_graph), alpha(_alpha), radius_of_steps(_radius_of_steps), bilateral(_bilateral) {
		infty = static_cast<cost_f_type>(this->graph.infty_value());
		list_of_prohib_moves;
	};

	virtual cost_f_type edge_cost(const int agent) const {
		return alpha * this->graph.strategy_of_player(agent).size();
	}

	virtual cost_f_type distance_cost(const int agent) const
	{ // Attention!!! If graph is disconnected, distance = infty
		auto	distances(this->graph.get_vector_of_dist(agent));

		if (!(this->graph.storage_dist()))
			distances = this->graph.vector_of_dist_weighted(agent);

		cost_f_type sum = 0.0;
		for (const auto& item : distances) {
			if (item != infty)
				sum += item;
			else
				return infty;
		}

		return static_cast<cost_f_type>(sum);
	}

	virtual cost_f_type cost(const int agent) const {
		return edge_cost(agent) + distance_cost(agent);
	}

	cost_f_type social_cost() const;

	virtual cost_f_type players_dist_cost_after_1_step(const int agent, const Step& step) const {
		auto distances = this->graph.distances_from_node_after_1_step(step, agent);
		cost_f_type sum = 0.0;
		for (const auto& item : distances) {
			if (item != infty)
				sum += item;
			else
				return infty;
		}
		return static_cast<cost_f_type>(sum);
	}

	virtual cost_f_type price_of_1_edge(const int v1, const int v2) const {
		if (this->graph.node_is_owned(v1, v2))
			return this->alpha;
		else
			return static_cast<cost_f_type>(0);
	}

	virtual cost_f_type players_edge_cost_after_1_step(const int agent, const Step& step) const {
		std::set<int> _curr_strategy = this->graph.strategy_of_player(agent);

		cost_f_type _edge_cost = edge_cost(agent);
		switch (step.step_name)
		{
		case 'a':
		{
			if (_curr_strategy.find(step.to) == _curr_strategy.end() && step.to != agent)
				_edge_cost += this->alpha;
			break;
		}
		case 'd':
		{
			if (_curr_strategy.find(step.to) == _curr_strategy.end() && step.from != agent)
				_edge_cost -= this->alpha;
			break;
		}
		default:
			break;
		}

		return  _edge_cost;
	}

	cost_f_type players_cost_after_1_step(const int agent, const Step& step) const {
		cost_f_type e_cost = players_edge_cost_after_1_step(agent, step);
		cost_f_type d_cost = players_dist_cost_after_1_step(agent, step);
		return e_cost + d_cost;
	}

	Step improving_response(const int agent, const bool random = 0) const;
		//greedy looking for improving response for agent, first check all additios, then deletions, then swaps
	// random means randomly chosen node to which strategy changes.

	Step pairwise_improving_response(const int agent, const bool _random) const;// to test
	//random = choose endpoint to which add or delete edge uar

	std::set<int> best_response(const int agent);

	Step best_greedy_response(const int agent) const;
	Step best_greedy_response_pairwise(const int agent) const;

	bool check_if_GE_unilateral(const bool output_impr_step = 0) const;

	//by playing imoroving response, agents are activated randomly (because it is faster)
	//outputs .dot-file with GE 
	//random parameter means how enpoints for new nodes are generated
	void greedy_find_GE(const std::string& extra_notes_in_output_filename, const bool round_robin, const bool random_endpoint, const bool _best_response, const int when_output_results, const bool output_graph_overview=1);
	void best_of_the_best_response_GE(const std::string& extra_notes_in_output_filename, const bool output_graph_overview, const int when_output_results);

	//agents are activated by round-robin, play improving responce
	bool find_IR_cycle(int limit_of_steps, bool output_grap_each_iter);

	//for tests
	void make_move(const int agent, const Step& step);

	void output_graph(std::string file_name) const;

	bool check_if_greedy_pairwise_stable(const bool output_impr_step = 0) const;

protected:
	Graph_type graph;
	cost_f_type alpha;
	cost_f_type beta = 0.0;
	int radius_of_steps;

	std::set<char> list_of_prohib_moves; // \subseteq of {'s','d','a'}

	bool bilateral;
	cost_f_type infty;

private:
	void increase_bitstring_by_1(std::string& bitstr, const std::set<int>& avoid_agent) {//avoids edges to itself
		for (auto i = 0; i < bitstr.size(); i++)
		{
			if (bitstr[i] == 1)
				bitstr[i] = 0;
			else if (avoid_agent.find(i) == avoid_agent.end())
			{
				bitstr[i] = 1;
				break;
			}
		}
	}
};



//==================  Homophilic NCG============


//this game is always bilateral
//we use friends-enemy concept, i.e., nodes are friends if they have a similar type (color)
class SchellingNCG : public Game<Node_colored_graph, double> {
private:
	int num_types = 2; 
	std::string timeline_file_name;
	std::vector<double> local_segregation;
public:
	explicit SchellingNCG(const Node_colored_graph& _input_graph, double _alpha, std::string _timeline_file_name) :
		Game<Node_colored_graph, double>(_input_graph, _alpha, 1), timeline_file_name(_timeline_file_name){
			local_segregation.resize(this->graph.graph_num_of_nodes(), 0.0);	
			for (int i = 0; i < this->graph.graph_num_of_nodes(); i++) {
				local_segregation[i] = static_cast<double>(static_cast<double>(this->graph.num_of_friends(i)) / this->graph.node_degree(i));
			}
	};

	explicit SchellingNCG(const Node_colored_graph& _input_graph, double _alpha, std::set<char>& _list_of_prohib_moves, std::string _timeline_file_name) :
		Game<Node_colored_graph, double>(_input_graph, _alpha, _list_of_prohib_moves, 1), timeline_file_name(_timeline_file_name) {
		local_segregation.resize(this->graph.graph_num_of_nodes(), 0.0);
		for (int i = 0; i < this->graph.graph_num_of_nodes(); i++) {
			local_segregation[i] = static_cast<double>(static_cast<double>(this->graph.num_of_friends(i)) / this->graph.node_degree(i));
		}
	};

	explicit SchellingNCG(const Node_colored_graph& _input_graph, double _alpha, int _radius_of_steps, std::set<char>& _list_of_prohib_moves, std::string _timeline_file_name) :
		Game<Node_colored_graph, double>(_input_graph, _alpha, _radius_of_steps, _list_of_prohib_moves, true), timeline_file_name(_timeline_file_name) {
		local_segregation.resize(this->graph.graph_num_of_nodes(), 0.0);
		for (int i = 0; i < this->graph.graph_num_of_nodes(); i++) {
			local_segregation[i] = static_cast<double>(static_cast<double>(this->graph.num_of_friends(i)) / this->graph.node_degree(i));
		}
	};

	double edge_cost(const int agent) const;
	double players_edge_cost_after_1_step(const int agent, const Step& step) const;
	double cost_of_the_neighborhood(const int _deg_value, const int _num_friends_value) const;

	Step best_greedy_response_pairwise(const int agent, const double approximation_factor) const;
	Step pairwise_improving_response(const int agent, const bool _random, const double approximation_factor) const;
	double greedy_find_GE(const std::string& extra_notes_in_output_filename, const bool round_robin, const bool random_endpoint, const bool _best_response, const int when_output_results, const bool output_graph_overview = 1, const double approximation_factor=1);

	double compute_avg_segr() const;
	void recompute_local_segr(int active_agent_1, int active_agent_2);
};


//============================================================
//****************** Implementation **************************
//============================================================
template <typename Graph_type, typename cost_f_type>
cost_f_type Game<Graph_type, cost_f_type>::social_cost() const {
	cost_f_type tot_cost = static_cast<cost_f_type>(0);
	for (auto i = 0; i < graph.graph_num_of_nodes(); i++)
		tot_cost += cost(i);
	return tot_cost;
}

template <typename Graph_type, typename cost_f_type>
Step Game<Graph_type, cost_f_type>::improving_response(const int agent, const bool random)  const //TODO randomization
//greedy looking for improving response for agent, first check all additions, then deletions, then swaps
//do not consider self-loops
{
	std::vector<int> accessible_nodes;
	std::set<int> k_neighbors;
	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		k_neighbors = graph.k_neighborhood(agent, radius_of_steps);
	}

	auto curr_dist_cost = distance_cost(agent); //need this for faster del-computations
	auto curr_agent_cost = edge_cost(agent) + curr_dist_cost;

	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		if (!graph.nodes_are_connected(agent, i) && i != agent) //for add-only version creating the list of accessible nodes takes too long !!!
			if (k_neighbors.size() == 0 || k_neighbors.find(i) != k_neighbors.end()) //if the game is non-local or the element is in k-neighborhood
				accessible_nodes.push_back(i);
	}

	//check additions
	if (list_of_prohib_moves.find('a') == list_of_prohib_moves.end()) {
		if (random) {
			std::set<int> list_of_still_accessible_elem(accessible_nodes.begin(), accessible_nodes.end());
			for (auto i = 0; i < accessible_nodes.size(); i++) {
				auto add_to = list_of_still_accessible_elem.begin();
				std::advance(add_to, rand() % list_of_still_accessible_elem.size()); //get k-th random element of the set 
				Step step('a', -1, *add_to);
				auto new_cost = players_cost_after_1_step(agent, step);
				if (new_cost < curr_agent_cost)
					return  step;
				list_of_still_accessible_elem.erase(add_to);
			}
		}
		else
			for (const auto& i : accessible_nodes) {
				Step step('a', -1, i);
				cost_f_type new_cost = players_cost_after_1_step(agent, step);
				if (new_cost < curr_agent_cost)
					return  step;
			}
	}

	//check deletions
	if (list_of_prohib_moves.find('d') == list_of_prohib_moves.end()) {
		for (const auto& i : graph.strategy_of_player(agent)) {
			Step step('d', i, -1);
			cost_f_type new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
			if (new_dist_cost < curr_dist_cost + price_of_1_edge(agent, i) && new_dist_cost != this->graph.infty_value()) {
				cost_f_type new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
				if (new_cost < curr_agent_cost)
					return step;
			}
		}
	}

	//check swaps
	if (list_of_prohib_moves.find('s') == list_of_prohib_moves.end()) {
		for (const auto& node_from : graph.strategy_of_player(agent)) {
			for (const auto& node_to : accessible_nodes) {
				Step step('s', node_from, node_to);
				cost_f_type new_cost = players_cost_after_1_step(agent, step);
				if (new_cost < curr_agent_cost)
					return step;
			}
		}
	}

	//if no improvement, return empty step
	return  Step(' ', -1, -1);	
}

template <typename Graph_type, typename cost_f_type>
Step Game<Graph_type, cost_f_type>::pairwise_improving_response(const int agent, const bool _random) const {	//add "random" parameter for the random activation scheme
	std::vector<int> accessible_nodes;
	std::set<int> k_neighborhood;
	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		k_neighborhood = graph.k_neighborhood(agent, radius_of_steps);
	}

	auto curr_dist_cost = distance_cost(agent); //need this for faster del-computations
	auto curr_agent_cost = edge_cost(agent) + curr_dist_cost;


	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		if (!graph.nodes_are_connected(agent, i) && i != agent)
			if (k_neighborhood.size() == 0 || std::binary_search(k_neighborhood.begin(), k_neighborhood.end(), i)) //if the game is local and the element is in k-neighborhood
				accessible_nodes.push_back(i);
	}

	if (_random) {
		std::set<int> list_of_not_checked_endpoints;

		if (radius_of_steps < graph.graph_num_of_nodes()) //if local
		{
			for (int i = 0; i < graph.graph_num_of_nodes(); i++)
				if (std::binary_search(k_neighborhood.begin(), k_neighborhood.end(), i) && i != agent)
					list_of_not_checked_endpoints.insert(i);
		}
		else
		{
			for (int i = 0; i < graph.graph_num_of_nodes(); i++)
			{
				if (i != agent)
					list_of_not_checked_endpoints.insert(i);
			}
		}

		while (list_of_not_checked_endpoints.size() > 0) {
			auto i = list_of_not_checked_endpoints.begin();
			std::advance(i, rand() % list_of_not_checked_endpoints.size()); // get i-th element of the set
			int active_endpoint = *i;

			list_of_not_checked_endpoints.erase(active_endpoint);

			if (this->list_of_prohib_moves.find('d') == this->list_of_prohib_moves.end()
				&& this->graph.nodes_are_connected(agent, active_endpoint))
			{//try deletion
				Step step('d', active_endpoint, -1);
				cost_f_type new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
				if (new_dist_cost < curr_dist_cost + price_of_1_edge(agent, active_endpoint) && new_dist_cost != this->graph.infty_value()) {
					cost_f_type new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
					if (new_cost < curr_agent_cost)
						return step;
				}
			}
			else if (this->list_of_prohib_moves.find('a') == this->list_of_prohib_moves.end())
			{//try to perform addition
				Step step('a', -1, active_endpoint);

				cost_f_type new_cost_agent_1 = players_cost_after_1_step(agent, step);
				cost_f_type cost_agent_2 = cost(active_endpoint);
				cost_f_type new_cost_agent_2 = players_cost_after_1_step(active_endpoint, Step('a', -1, agent));
				if (new_cost_agent_1 < curr_agent_cost && new_cost_agent_2 < cost_agent_2)
					return  step;
			}
		}
	}
	else {
		//check additions
		if (list_of_prohib_moves.find('a') == list_of_prohib_moves.end()) {
			for (const auto& i : accessible_nodes) {
				Step step('a', -1, i);

				cost_f_type new_cost_agent_1 = players_cost_after_1_step(agent, step);
				cost_f_type cost_agent_2 = cost(i);
				cost_f_type new_cost_agent_2 = players_cost_after_1_step(i, Step('a', -1, agent));
				if (new_cost_agent_1 < curr_agent_cost && new_cost_agent_2 < cost_agent_2)
					return  step;
			}
		}

		//check deletions
		if (list_of_prohib_moves.find('d') == list_of_prohib_moves.end()) {
			for (const auto& i : graph.strategy_of_player(agent)) {
				Step step('d', i, -1);
				cost_f_type new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
				if (new_dist_cost < curr_dist_cost + price_of_1_edge(agent, i) && new_dist_cost != this->graph.infty_value()) {
					cost_f_type new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
					if (new_cost < curr_agent_cost)
						return step;
				}
			}
		}
	}

	//if no improvement, return empty step
	return  Step(' ', -1, -1);
}



template <typename Graph_type, typename cost_f_type>
std::set<int> Game<Graph_type, cost_f_type>::best_response(const int agent) {		//better variant???
	//greedy check all possible strategies
	std::set<int> best_strategy = {};
	graph.set_strategy(agent, best_strategy);
	cost_f_type best_cost = cost(agent);

	std::string bitstring(graph.graph_num_of_nodes(), 0);
	std::string final_str(graph.graph_num_of_nodes(), 1);

	std::set<int> restrict_set = graph.neighborhood(agent); //do not consider edges to this nodes because they are already in the graph
	restrict_set.insert(agent);

	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		std::set<int> k_neighborhood;
		k_neighborhood = graph.k_neighborhood(agent, radius_of_steps);
		for (int i = 0; i < graph.graph_num_of_nodes(); i++)
			if (k_neighborhood.find(i) != k_neighborhood.end())
				restrict_set.insert(i);
	}

	for (auto& i : restrict_set)
		final_str[i] = 0;

	do {
		increase_bitstring_by_1(bitstring, restrict_set);

		std::set<int> new_strategy;
		for (auto i = 0; i < bitstring.size(); i++)
			if (bitstring[i])
				new_strategy.insert(i);

		graph.set_strategy(agent, new_strategy);

		cost_f_type new_cost = cost(agent);
		if (new_cost < best_cost) {
			best_cost = new_cost;
			best_strategy = new_strategy;
		}
	} while (bitstring != final_str);

	return best_strategy;
}


template <typename Graph_type, typename cost_f_type>
Step Game<Graph_type, cost_f_type>::best_greedy_response(const int agent) const {
	std::vector<int> accessible_nodes;
	std::set<int> k_neighborhood;
	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		k_neighborhood = graph.k_neighborhood(agent, radius_of_steps);
	}

	cost_f_type curr_agent_cost = cost(agent);

	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		if (graph.neighborhood(agent).count(i) == 0 && i != agent)
			if (k_neighborhood.size() == 0 || std::binary_search(k_neighborhood.begin(), k_neighborhood.end(), i)) //if the game is local and the element is in k-neighborhood
				accessible_nodes.push_back(i);
	}

	Step best_move(' ', -1, -1);

	if (this->list_of_prohib_moves.find('a') == list_of_prohib_moves.end()) {
		//check additions
		for (const auto& i : accessible_nodes) {
			cost_f_type new_cost = players_cost_after_1_step(agent, Step('a', -1, i));

			if (new_cost < curr_agent_cost)
			{
				best_move = Step('a', -1, i);
				curr_agent_cost = new_cost;
				std::cout << "add " << agent << ' ' << i << ", dist=" << distance_cost(agent) << ", edge_cost=" << edge_cost(agent) << ", new cost = " << new_cost << '\n';
			}
		}
	}

	if (this->list_of_prohib_moves.find('d') == this->list_of_prohib_moves.end()) {
		//check deletions
		for (const auto& i : graph.strategy_of_player(agent)) {
			cost_f_type new_cost = players_cost_after_1_step(agent, Step('d', i, -1));
			if (new_cost < curr_agent_cost)
			{
				best_move = Step('d', i, -1);
				curr_agent_cost = new_cost;
				std::cout << "del " << agent << ' ' << i << ", dist=" << distance_cost(agent) << ", edge_cost=" << edge_cost(agent) << ", new cost = " << new_cost << '\n';
			}
		}
	}

	if (this->list_of_prohib_moves.find('s') == this->list_of_prohib_moves.end()) {
		//check swaps
		for (const auto& node_from : graph.strategy_of_player(agent)) {
			for (const auto& node_to : accessible_nodes) {
				cost_f_type new_cost = players_cost_after_1_step(agent, Step('s', node_from, node_to));
				if (new_cost < curr_agent_cost)
				{
					best_move = Step('s', node_from, node_to);
					curr_agent_cost = new_cost;
					std::cout << "swap " << agent << ' ' << node_from << " to: " << node_to << ", dist=" << distance_cost(agent) << ", edge_cost=" << edge_cost(agent) << ", new cost = " << new_cost << '\n';
				}
			}
		}
	}

	//if no improvement, return empty step
	return  best_move;														

}

template<typename Graph_type, typename cost_f_type>
Step Game<Graph_type, cost_f_type>::best_greedy_response_pairwise(const int agent) const {
	std::vector<int> accessible_nodes;
	std::set<int> k_neighborhood;
	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		k_neighborhood = graph.k_neighborhood(agent, radius_of_steps);
	}

	cost_f_type curr_dist_cost = distance_cost(agent);
	cost_f_type curr_agent_cost = curr_dist_cost + edge_cost(agent);

	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		if (graph.neighborhood(agent).count(i) == 0 && i != agent)
			if (k_neighborhood.size() == 0 || std::binary_search(k_neighborhood.begin(), k_neighborhood.end(), i)) //if the game is local and the element is in k-neighborhood
				accessible_nodes.push_back(i);
	}

	Step best_move(' ', -1, -1);

	if (this->list_of_prohib_moves.find('a') == this->list_of_prohib_moves.end()) {
		//check additions
		for (const auto& i : accessible_nodes) {
			Step step('a', -1, i);
			cost_f_type new_dist_cost_agent_1 = players_dist_cost_after_1_step(agent, step);
			cost_f_type new_cost_agent_1 = new_dist_cost_agent_1 + players_edge_cost_after_1_step(agent, step);
			cost_f_type cost_agent_2 = cost(i);
			cost_f_type new_cost_agent_2 = players_cost_after_1_step(i, Step('a', -1, agent));

			if (new_cost_agent_1 < curr_agent_cost && new_cost_agent_2 < cost_agent_2)

			{
				best_move = Step('a', -1, i);
				curr_agent_cost = new_cost_agent_1;
			}
		}
	}

	if (this->list_of_prohib_moves.find('d') == this->list_of_prohib_moves.end()) {
		//check deletions
		for (const auto& active_endpoint : graph.strategy_of_player(agent)) {
			Step step('d', active_endpoint, -1);
			cost_f_type new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
			if (new_dist_cost < curr_dist_cost + price_of_1_edge(agent, active_endpoint) && new_dist_cost != this->graph.infty_value()) {
				cost_f_type new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
				if (new_cost < curr_agent_cost) {
					best_move = Step('d', active_endpoint, -1);
					curr_agent_cost = new_cost;
					curr_dist_cost = new_dist_cost;
					std::cout << "del " << agent << ' ' << active_endpoint << '\n';
				}
			}
		}
	}

	//if no improvement, return empty step
	return  best_move;														
}

template <typename Graph_type, typename cost_f_type>
bool Game<Graph_type, cost_f_type>::check_if_GE_unilateral(const bool output_impr_step) const {
	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		{
			Step move;
			move = improving_response(i);

			if (move.step_name != ' ')
			{
				if (output_impr_step) std::cout << i << " has improving move: " << move;
				return 0;
			}
		}
	}
	return 1;
}

template <typename Graph_type, typename cost_f_type>
void Game<Graph_type, cost_f_type>::greedy_find_GE(const std::string& extra_notes_in_output_filename, const bool round_robin, const bool random_endpoint, const bool _best_response, const int when_output_results, const bool output_graph_overview) {

	std::set<int> set_of_nonactive_agents, full_set;
	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		set_of_nonactive_agents.insert(i);
		full_set.insert(i);
	}
	int count = 0;
	int active_agent = -1;

	while (set_of_nonactive_agents.size() != 0) {
		auto i = set_of_nonactive_agents.begin();

		if (round_robin)
			active_agent = (active_agent + 1) % (this->graph.graph_num_of_nodes());
		else
		{
			std::advance(i, rand() % set_of_nonactive_agents.size()); // get i-th random element of the set
			active_agent = *i;
		}

		Step move;
		if (this->bilateral)
			if (_best_response)
				move = best_greedy_response_pairwise(active_agent);
			else
				move = pairwise_improving_response(active_agent, random_endpoint);
		else
			if (_best_response)
				move = best_greedy_response(active_agent);
			else
				move = improving_response(active_agent, random_endpoint);

		if (move.step_name == ' ')
		{
			set_of_nonactive_agents.erase(active_agent);
			std::cout << "agent " << active_agent << " is happy\n";
		}
		else {
			count++;
			this->graph.perform_step(active_agent, move, this->bilateral, 1);

			if ((when_output_results > 0) && (count % when_output_results) == 0) {
				//output degree distribution
				std::vector<int> deg_distribution(this->graph.graph_num_of_nodes(), 0);
				for (int i = 0; i < this->graph.graph_num_of_nodes(); i++)
					deg_distribution[this->graph.node_degree(i)]++;

				std::ofstream out_deg_distr_file("alpha_" + std::to_string(this->alpha) + "_step_" + std::to_string(count) + extra_notes_in_output_filename + "_deg_distribution.dot");
				if (out_deg_distr_file.is_open())
				{
					for (const auto i : deg_distribution)
						out_deg_distr_file << std::to_string(i) << '\n';
					out_deg_distr_file.close();
				}
				else
				{
					std::cout << "Cannot open out_file for the degree_distribution output";
					exit(0);
				}

				graph.output_graph_to_dot_file("alpha_" + std::to_string(this->alpha) + "_step_" + std::to_string(count) + extra_notes_in_output_filename + ".dot");
			}
			set_of_nonactive_agents = full_set;					//check this for better variant
			std::cout << "agent " << active_agent << " does " << move << '\n';
		}
	}
	cost_f_type cost = social_cost();

	std::string extra_notes = "";
	if (radius_of_steps < graph.graph_num_of_nodes()) {
		extra_notes += "_" + std::to_string(radius_of_steps) + "_local";
	}
	else
		extra_notes += "_global";

	extra_notes += "_" + std::to_string(graph.avg_clustering()) + "_avg_clust_" + std::to_string(this->graph.diameter()) + "_D";
	this->graph.output_graph_to_dot_file("GE_" + std::to_string(graph.graph_num_of_nodes()) + "_n_" + extra_notes_in_output_filename + extra_notes + ".dot");

	//output all statistics for the resulting graph
	if (output_graph_overview) {
		//1. degree distribution
		std::vector<int> deg_distribution(this->graph.graph_num_of_nodes(), 0);
		for (int i = 0; i < this->graph.graph_num_of_nodes(); i++)
			deg_distribution[this->graph.node_degree(i)]++;

		std::ofstream out_deg_distr_file(std::to_string(graph.graph_num_of_nodes()) + "_n_" + extra_notes_in_output_filename + "_deg_distribution.dot");
		if (out_deg_distr_file.is_open())
		{
			for (const auto i : deg_distribution)
				out_deg_distr_file << std::to_string(i) << '\n';
			out_deg_distr_file.close();
		}
		else
		{
			std::cout << "Cannot open out_file for the degree_distribution output";
			exit(0);
		}


		//3. degree x local CC
		std::ofstream out_deg_CC_file(std::to_string(graph.graph_num_of_nodes()) + "_n_" + extra_notes_in_output_filename + "_node_deg_local_CC.dot");
		if (out_deg_CC_file.is_open())
		{
			for (int i = 0; i < this->graph.graph_num_of_nodes(); i++) {

				out_deg_CC_file << this->graph.node_degree(i) << '	' << this->graph.local_clustering(i) << '\n';
			}
			out_deg_CC_file.close();
		}
		else
		{
			std::cout << "Cannot open out_file for the deg x CC output";
			exit(0);
		}
	}
}

template <typename Graph_type, typename cost_f_type>
void Game<Graph_type, cost_f_type>::best_of_the_best_response_GE(const std::string& extra_notes_in_output_filename, const bool output_graph_overview, const int when_output_results) {

	int count = 0;
	Step best_BR;
	int BR_active_agent;
	cost_f_type best_delta_cost;

	do {
		std::vector<std::pair<cost_f_type, cost_f_type>> list_of_agents_costs;
		for (int agent = 0; agent < this->graph.graph_num_of_nodes(); agent++) {
			list_of_agents_costs.push_back(std::pair<cost_f_type, cost_f_type>(edge_cost(agent), distance_cost(agent)));
		}

		best_delta_cost = 0.0;

		for (int active_agent_1 = 0; active_agent_1 < this->graph.graph_num_of_nodes(); active_agent_1++) {
			if (this->bilateral)
				for (int active_agent_2 = active_agent_1 + 1; active_agent_2 < this->graph.graph_num_of_nodes(); active_agent_2++)
				{
					if (this->graph.nodes_are_connected(active_agent_1, active_agent_2) && (this->list_of_prohib_moves.find('d') == this->list_of_prohib_moves.end())) //then try delete the edge by both endpoints sequentialy
					{
						Step move('d', active_agent_2, -1);
						cost_f_type new_dist_cost = this->players_dist_cost_after_1_step(active_agent_1, move);
						if (new_dist_cost < list_of_agents_costs[active_agent_1].second + price_of_1_edge(active_agent_1, active_agent_2) && new_dist_cost != this->graph.infty_value()) {
							cost_f_type new_delta_cost = players_edge_cost_after_1_step(active_agent_1, move) + new_dist_cost - (list_of_agents_costs[active_agent_1].first + list_of_agents_costs[active_agent_1].second);

							if (new_delta_cost < best_delta_cost) {
								best_BR = move;
								BR_active_agent = active_agent_1;
								best_delta_cost = new_delta_cost;
								std::cout << "new BR: del " << std::to_string(active_agent_1) << "->" << std::to_string(active_agent_2) << ", new delta_cost = " << best_delta_cost << '\n';
							}
						}
						//check del by the second endpoint
						move = Step('d', active_agent_1, -1);
						new_dist_cost = this->players_dist_cost_after_1_step(active_agent_2, move);
						if (new_dist_cost < list_of_agents_costs[active_agent_2].second + price_of_1_edge(active_agent_2, active_agent_1) && new_dist_cost != this->graph.infty_value()) {
							cost_f_type new_delta_cost = players_edge_cost_after_1_step(active_agent_2, move) + new_dist_cost - (list_of_agents_costs[active_agent_2].first + list_of_agents_costs[active_agent_2].second);

							if (new_delta_cost < best_delta_cost) {
								best_BR = move;
								BR_active_agent = active_agent_2;
								best_delta_cost = new_delta_cost;
								std::cout << "new BR: del " << std::to_string(active_agent_2) << "->" << std::to_string(active_agent_1) << ", new delta_cost = " << best_delta_cost << '\n';
							}
						}
					}
					else if (!(this->graph.nodes_are_connected(active_agent_1, active_agent_2)) && (this->list_of_prohib_moves.find('a') == this->list_of_prohib_moves.end())) //try addition of the missing edge
					{
						Step move('a', -1, active_agent_2);
						cost_f_type new_delta_cost_agent_1 = players_cost_after_1_step(active_agent_1, move) - (list_of_agents_costs[active_agent_1].first + list_of_agents_costs[active_agent_1].second);
						cost_f_type new_delta_cost_agent_2 = players_cost_after_1_step(active_agent_2, Step('a', -1, active_agent_1)) - (list_of_agents_costs[active_agent_2].first + list_of_agents_costs[active_agent_2].second);
						if (new_delta_cost_agent_1 < 0 && new_delta_cost_agent_2 < 0 && std::min(new_delta_cost_agent_1, new_delta_cost_agent_2) < best_delta_cost)
						{
							best_BR = move;
							BR_active_agent = active_agent_1; //since bilateral version, we don't care who is the owner of the best edge
							best_delta_cost = std::min(new_delta_cost_agent_1, new_delta_cost_agent_2);
							std::cout << "new BR: add " << std::to_string(active_agent_1) << "->" << std::to_string(active_agent_2) << ", new delta_cost = " << best_delta_cost << '\n';
						}

					}
				}
			else
				for (int active_agent_2 = 0; active_agent_2 < this->graph.graph_num_of_nodes(); active_agent_2++)
				{
					if (this->graph.nodes_are_connected(active_agent_1, active_agent_2) && (this->list_of_prohib_moves.find('d') == this->list_of_prohib_moves.end())) //then try delete the edge by both endpoints sequentialy
					{
						Step move('d', active_agent_2, -1);
						cost_f_type new_dist_cost = this->players_dist_cost_after_1_step(active_agent_1, move);
						if (new_dist_cost < list_of_agents_costs[active_agent_1].second + price_of_1_edge(active_agent_1, active_agent_2) && new_dist_cost != this->graph.infty_value()) {
							cost_f_type new_delta_cost = players_edge_cost_after_1_step(active_agent_1, move) + new_dist_cost - (list_of_agents_costs[active_agent_1].first + list_of_agents_costs[active_agent_1].second);
							if (new_delta_cost < best_delta_cost) {
								best_BR = move;
								BR_active_agent = active_agent_1;
								best_delta_cost = new_delta_cost;
								std::cout << "new BR: del " << std::to_string(active_agent_1) << "->" << std::to_string(active_agent_2) << ", new delta_cost = " << best_delta_cost << '\n';
							}
						}
						//check del by the second endpoint
						move = Step('d', active_agent_1, -1);
						new_dist_cost = this->players_dist_cost_after_1_step(active_agent_2, move);
						if (new_dist_cost < list_of_agents_costs[active_agent_2].second + price_of_1_edge(active_agent_2, active_agent_1) && new_dist_cost != this->graph.infty_value()) {
							cost_f_type new_delta_cost = players_edge_cost_after_1_step(active_agent_2, move) + new_dist_cost - (list_of_agents_costs[active_agent_2].first + list_of_agents_costs[active_agent_2].second);
							if (new_delta_cost < best_delta_cost) {
								best_BR = move;
								BR_active_agent = active_agent_2;
								best_delta_cost = new_delta_cost;
								std::cout << "new BR: del " << std::to_string(active_agent_2) << "->" << std::to_string(active_agent_1) << ", new delta_cost = " << best_delta_cost << '\n';
							}
						}
					}
					else if (!(this->graph.nodes_are_connected(active_agent_1, active_agent_2)) && (this->list_of_prohib_moves.find('a') == this->list_of_prohib_moves.end()))//try addition of the missing edge
					{
						Step move('a', -1, active_agent_2);
						cost_f_type new_delta_cost_agent_1 = players_cost_after_1_step(active_agent_1, move) - (list_of_agents_costs[active_agent_1].first + list_of_agents_costs[active_agent_1].second);

						if (new_delta_cost_agent_1 < best_delta_cost)
						{
							best_BR = move;
							BR_active_agent = active_agent_1;
							best_delta_cost = new_delta_cost_agent_1;
							std::cout << "new BR: add " << std::to_string(active_agent_1) << "->" << std::to_string(active_agent_2) << ", new delta_cost = " << best_delta_cost << '\n';
						}

						move = Step('a', -1, active_agent_1);
						cost_f_type new_delta_cost_agent_2 = players_cost_after_1_step(active_agent_2, move) - (list_of_agents_costs[active_agent_2].first + list_of_agents_costs[active_agent_2].second);
						if (new_delta_cost_agent_2 < best_delta_cost)
						{
							best_BR = move;
							BR_active_agent = active_agent_2;
							best_delta_cost = new_delta_cost_agent_2;
							std::cout << "new BR: add " << std::to_string(active_agent_2) << "->" << std::to_string(active_agent_2) << ", new delta_cost = " << best_delta_cost << '\n';
						}

					}
				}
		}

		if (best_delta_cost != 0) {
			count++;
			this->graph.perform_step(BR_active_agent, best_BR, this->bilateral, 1);
			std::cout << "WINNER: " << std::to_string(BR_active_agent) << ", " << best_BR << '\n';

			if ((when_output_results > 0) && (count % when_output_results) == 0) {

				//output degree distribution
				std::vector<int> deg_distribution(this->graph.graph_num_of_nodes(), 0);
				for (int i = 0; i < this->graph.graph_num_of_nodes(); i++)
					deg_distribution[this->graph.node_degree(i)]++;

				std::ofstream out_deg_distr_file("alpha_" + std::to_string(this->alpha) + "_step_" + std::to_string(count) + extra_notes_in_output_filename + "_deg_distribution.dot");
				if (out_deg_distr_file.is_open())
				{
					for (const auto i : deg_distribution)
						out_deg_distr_file << std::to_string(i) << '\n';
					out_deg_distr_file.close();
				}
				else
				{
					std::cout << "Cannot open out_file for the degree_distribution output";
					exit(0);
				}

				graph.output_graph_to_dot_file("alpha_" + std::to_string(this->alpha) + "_step_" + std::to_string(count) + extra_notes_in_output_filename + ".dot");
			}
		}
	} while (best_delta_cost != 0);


	cost_f_type cost = social_cost();

	std::string extra_notes = "";
	if (radius_of_steps < graph.graph_num_of_nodes()) {
		extra_notes += "_" + std::to_string(radius_of_steps) + "_local";
	}
	else
		extra_notes += "_global";

	extra_notes += "_" + std::to_string(graph.avg_clustering()) + "_avg_clust_" + std::to_string(this->graph.diameter()) + "_D";
	this->graph.output_graph_to_dot_file("GE_" + std::to_string(graph.graph_num_of_nodes()) + "_n_" + extra_notes_in_output_filename + extra_notes + ".dot");

	//output all statistics for the resulting graph

	if (output_graph_overview) {
		//1. degree distribution
		std::vector<int> deg_distribution(this->graph.graph_num_of_nodes(), 0);
		for (int i = 0; i < this->graph.graph_num_of_nodes(); i++)
			deg_distribution[this->graph.node_degree(i)]++;

		std::ofstream out_deg_distr_file(std::to_string(graph.graph_num_of_nodes()) + "_n_" + extra_notes_in_output_filename + "_deg_distribution.dot");
		if (out_deg_distr_file.is_open())
		{
			for (const auto i : deg_distribution)
				out_deg_distr_file << std::to_string(i) << '\n';
			out_deg_distr_file.close();
		}
		else
		{
			std::cout << "Cannot open out_file for the degree_distribution output";
			exit(0);
		}


		//3. degree x local CC
		std::ofstream out_deg_CC_file(std::to_string(graph.graph_num_of_nodes()) + "_n_" + extra_notes_in_output_filename + "_node_deg_local_CC.dot");
		if (out_deg_CC_file.is_open())
		{
			for (int i = 0; i < this->graph.graph_num_of_nodes(); i++) {

				out_deg_CC_file << this->graph.node_degree(i) << '	' << this->graph.local_clustering(i) << '\n';
			}
			out_deg_CC_file.close();
		}
		else
		{
			std::cout << "Cannot open out_file for the deg x CC output";
			exit(0);
		}
	}
}

template <typename Graph_type, typename cost_f_type>
bool Game<Graph_type, cost_f_type>::find_IR_cycle(int limit_of_steps, bool output_grap_each_iter) {
	int iteration = 0;
	int active_agent = 0;

	std::vector<cost_f_type> list_of_cost_values;
	std::vector<cost_f_type> tmp_list;
	while (active_agent < graph.graph_num_of_nodes() && iteration < limit_of_steps) {
		Step impr_move;
		if (this->bilateral) impr_move = pairwise_improving_response(active_agent, 0);
		else impr_move = improving_response(active_agent, 0);

		if (impr_move.step_name != ' ')//if there is an improvement - improve
		{
			graph.perform_step(active_agent, impr_move, this->bilateral, 1);
			active_agent = 0;
			iteration++;
			cost_f_type upd_cost = social_cost();
			if (output_grap_each_iter)
				graph.output_graph_to_dot_file(std::to_string(iteration) + "_iter" + std::to_string(upd_cost) + ".dot");

			auto find_result = std::find(list_of_cost_values.begin(), list_of_cost_values.end(), upd_cost);
			if (find_result != list_of_cost_values.end()) {
				list_of_cost_values.push_back(upd_cost);
			}
			else {
				if (tmp_list.size() == 0)
					tmp_list.push_back(upd_cost);
				else {
					if (tmp_list[0] == upd_cost) {
						this->find_IR_cycle(tmp_list.size(), true); //thus, there is a cycle. So we restart our search from the current graph and do only num_of_iter==size_of_the_cycle
						return true;
					}
					else
					{
						tmp_list.push_back(upd_cost);
						if (std::search(list_of_cost_values.begin(), list_of_cost_values.end(), tmp_list.begin(), tmp_list.end()) != list_of_cost_values.end()) {
							list_of_cost_values.insert(list_of_cost_values.end(), tmp_list.begin(), tmp_list.end());
							tmp_list.clear();
						}
					}
				}
			}
		}
		else
			active_agent++;
	}
	return false;
}


template <typename Graph_type, typename cost_f_type>
void Game<Graph_type, cost_f_type>::make_move(const int agent, const Step& step) {
	graph.perform_step(agent, step, this->bilateral, 1);
}


template <typename Graph_type, typename cost_f_type>
void Game<Graph_type, cost_f_type>::output_graph(std::string file_name) const {
	graph.output_graph_to_dot_file(file_name);
}

template <typename Graph_type, typename cost_f_type>
bool Game<Graph_type, cost_f_type>::check_if_greedy_pairwise_stable(const bool output_impr_step) const {
	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		{
			Step move = pairwise_improving_response(i, 0);

			if (move.step_name != ' ')
			{
				if (output_impr_step) std::cout << i << " has improving move: " << move;
				return 0;
			}
		}
	}
	return 1;
}