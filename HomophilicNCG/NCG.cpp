#include "NCG.h"

double harmonic_number(int n) {
	double result = 0.0;
	for (int i = 1; i <= n; i++)
		result += (1.0 / i);
	return result;
}

double limit_comp_error(double val) {
	if (abs(val) < 0.00000001)
		return 0.0;
	else 
		return val;
}


double SchellingNCG::cost_of_the_neighborhood(const int _deg_value, const int _num_friends_value) const {
			//Harmonic enemies
	//return this->alpha * _deg_value + this->alpha * harmonic_number(_deg_value - _num_friends_value);

	//initial game
	return this->alpha * (double)_deg_value * (1.0 + 1.0/((double)_num_friends_value+1.0));
}


double SchellingNCG::edge_cost(const int agent) const {
	return this->cost_of_the_neighborhood(this->graph.node_degree(agent), this->graph.num_of_friends(agent));
}

double SchellingNCG::players_edge_cost_after_1_step(const int agent, const Step& step) const {
	std::set<int> _curr_strategy = this->graph.strategy_of_player(agent);

	double _edge_cost = 0.0;
	switch (step.step_name)
	{
	case 'a':
	{
		if (_curr_strategy.find(step.to) == _curr_strategy.end() && step.to != agent) {
			if (this->graph.agents_are_friends(agent, step.to))
				_edge_cost = this->cost_of_the_neighborhood(this->graph.node_degree(agent) + 1, this->graph.num_of_friends(agent) + 1); //this->alpha * (this->graph.node_degree(agent)+1.0) * (1.0/ (this->graph.num_of_friends(agent)+1.0) + 1.0);
			else
				_edge_cost = this->cost_of_the_neighborhood(this->graph.node_degree(agent) + 1, this->graph.num_of_friends(agent));
		}
		break;
	}
	case 'd':
	{
		if (_curr_strategy.find(step.to) == _curr_strategy.end() && step.from != agent)

			if (this->graph.agents_are_friends(agent, step.from))
				_edge_cost = this->cost_of_the_neighborhood(this->graph.node_degree(agent) - 1, this->graph.num_of_friends(agent) - 1);
			else
				_edge_cost = this->cost_of_the_neighborhood(this->graph.node_degree(agent) - 1, this->graph.num_of_friends(agent)); //_edge_cost -= this->alpha * (this->graph.num_of_friends(agent) - 1.0);
		break;
	}
	default:
		break;
	}

	return  _edge_cost;
}

Step SchellingNCG::best_greedy_response_pairwise(const int agent, const double approximation_factor) const {
	std::vector<int> accessible_nodes;
	std::set<int> k_neighborhood;
	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		k_neighborhood = graph.k_neighborhood(agent, radius_of_steps);
	}

	double curr_dist_cost = distance_cost(agent);
	double curr_agent_cost = curr_dist_cost + edge_cost(agent);

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
			double new_dist_cost_agent_1 = players_dist_cost_after_1_step(agent, step);
			double new_cost_agent_1 = new_dist_cost_agent_1 + players_edge_cost_after_1_step(agent, step);
			double cost_agent_2 = cost(i);
			double new_cost_agent_2 = players_cost_after_1_step(i, Step('a', -1, agent));

			if (limit_comp_error(new_cost_agent_1 * approximation_factor - curr_agent_cost) < 0.0 && limit_comp_error(new_cost_agent_2 * approximation_factor - cost_agent_2) < 0.0)

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
			double new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
			if (new_dist_cost != this->graph.infty_value()) {
				double new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
				if (limit_comp_error(new_cost * approximation_factor - curr_agent_cost) < 0.0) {
					best_move = Step('d', active_endpoint, -1);
					curr_agent_cost = new_cost;
					curr_dist_cost = new_dist_cost;
					//std::cout << "del " << agent << ' ' << active_endpoint << ", dist=" << curr_dist_cost << ", edge_cost=" << edge_cost(agent) << ", new cost = " << new_cost << '\n';
					std::cout << "del " << agent << ' ' << active_endpoint << '\n';
				}
			}
		}
	}


	//if no improvement, return empty step
	return  best_move;
}

Step SchellingNCG::pairwise_improving_response(const int agent, const bool _random, const double approximation_factor) const {
	std::vector<int> accessible_nodes;
	std::set<int> k_neighborhood;
	if (radius_of_steps < graph.graph_num_of_nodes()) { //if the game is local
		k_neighborhood = graph.k_neighborhood(agent, radius_of_steps);
	}

	auto curr_dist_cost = distance_cost(agent); //need this for faster del-computations
	auto curr_agent_cost = edge_cost(agent) + curr_dist_cost;


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
			bool add_only = 0;
			if (this->list_of_prohib_moves.find('d') != this->list_of_prohib_moves.end()) //if add-only
			{
				add_only = 1;
			}

			for (int i = 0; i < graph.graph_num_of_nodes(); i++)
			{
				if (i != agent && (!add_only || !this->graph.nodes_are_connected(agent, i)))
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
				double new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
				if (new_dist_cost != this->graph.infty_value()) {
					double new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
					if (limit_comp_error(new_cost * approximation_factor - curr_agent_cost) < 0.0)
						return step;
				}
			}
			else if (this->list_of_prohib_moves.find('a') == this->list_of_prohib_moves.end()
				&& !this->graph.nodes_are_connected(agent, active_endpoint))
			{//try to perform addition
				Step step('a', -1, active_endpoint);

				double new_cost_agent_1 = players_cost_after_1_step(agent, step);
				double cost_agent_2 = cost(active_endpoint);
				double new_cost_agent_2 = players_cost_after_1_step(active_endpoint, Step('a', -1, agent));
				if (limit_comp_error(new_cost_agent_1 * approximation_factor - curr_agent_cost) < 0.0 && limit_comp_error(new_cost_agent_2 * approximation_factor - cost_agent_2) < 0.0)
					return  step;
			}
		}
	}
	else {								//if we need to check every step in a fixed order

		for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
			if (!graph.nodes_are_connected(agent, i) && i != agent)
				if (k_neighborhood.size() == 0 || std::binary_search(k_neighborhood.begin(), k_neighborhood.end(), i)) //if the game is local and the element is in k-neighborhood
					accessible_nodes.push_back(i);
		}

		//check additions
		if (list_of_prohib_moves.find('a') == list_of_prohib_moves.end()) {
			for (const auto& i : accessible_nodes) {
				Step step('a', -1, i);

				double  new_cost_agent_1 = players_cost_after_1_step(agent, step);
				double  cost_agent_2 = cost(i);
				double  new_cost_agent_2 = players_cost_after_1_step(i, Step('a', -1, agent));
				if (limit_comp_error(new_cost_agent_1 * approximation_factor - curr_agent_cost) < 0.0 && limit_comp_error(new_cost_agent_2 * approximation_factor - cost_agent_2) < 0.0)
					return  step;
			}
		}

		//check deletions
		if (list_of_prohib_moves.find('d') == list_of_prohib_moves.end()) {
			for (const auto& i : graph.strategy_of_player(agent)) {
				Step step('d', i, -1);
				double  new_dist_cost = this->players_dist_cost_after_1_step(agent, step);
				if (new_dist_cost != this->graph.infty_value()) {
					double  new_cost = players_edge_cost_after_1_step(agent, step) + new_dist_cost;
					if (limit_comp_error(new_cost * approximation_factor - curr_agent_cost) < 0.0)
						return step;
				}
			}
		}
	}

	//if no improvement, return empty step
	return  Step(' ', -1, -1);
}



double  SchellingNCG::compute_avg_segr() const {
	double result = 0.0;
	for (int i = 0; i < this->graph.graph_num_of_nodes(); i++)
		result += this->local_segregation[i];

	return result / this->graph.graph_num_of_nodes();
}

void SchellingNCG::recompute_local_segr(int active_agent_1, int active_agent_2)
{
	local_segregation[active_agent_1] = static_cast<double>(static_cast<double>(this->graph.num_of_friends(active_agent_1)) / this->graph.node_degree(active_agent_1));
	local_segregation[active_agent_2] = static_cast<double>(static_cast<double>(this->graph.num_of_friends(active_agent_2)) / this->graph.node_degree(active_agent_2));
}

double SchellingNCG::greedy_find_GE(const std::string& extra_notes_in_output_filename, const bool round_robin, const bool random_endpoint, const bool _best_response, const int when_output_results, const bool output_graph_overview, const double approximation_factor)
{
	std::set<int> set_of_nonactive_agents, full_set;
	for (auto i = 0; i < graph.graph_num_of_nodes(); i++) {
		set_of_nonactive_agents.insert(i);
		full_set.insert(i);
	}
	int count = 0;//number of steps

	int active_agent = -1;

	double avg_segr = compute_avg_segr();


	std::ofstream timeline_avg_segr_file(this->timeline_file_name + "_avg_segr_timeline.txt", std::ios::out | std::ios::app);
	std::ofstream timeline_moves_file(this->timeline_file_name + "_moves_timeline.txt", std::ios::out | std::ios::app);
	int current_num_of_edges = this->graph.num_of_edges();
	timeline_moves_file << current_num_of_edges << "	initial state";
	timeline_moves_file.close();

	if (timeline_avg_segr_file.is_open())
		timeline_avg_segr_file << avg_segr;
	else {
		std::cout << "Cannot open timeline_avg_segr_file for the time line output";
		exit(0);
	}
	timeline_avg_segr_file.close();

	while (set_of_nonactive_agents.size() != 0 && count <= 50'000) {
		auto i = set_of_nonactive_agents.begin();

		if (round_robin)
			active_agent = (active_agent + 1) % (this->graph.graph_num_of_nodes());
		else
		{
			std::advance(i, rand() % set_of_nonactive_agents.size()); // get i-th random element of the set
			active_agent = *i;
		}

		Step move;

		if (_best_response)
			move = best_greedy_response_pairwise(active_agent, approximation_factor);
		else
			move = pairwise_improving_response(active_agent, random_endpoint, approximation_factor);


		if (move.step_name == ' ')
		{
			set_of_nonactive_agents.erase(active_agent);
			std::cout << "agent " << active_agent << " is happy\n";
		}
		else {
			count++;
			this->graph.perform_step(active_agent, move, this->bilateral, 1);

			set_of_nonactive_agents = full_set;					//check this for better variant
			std::cout << "agent " << active_agent << " does " << move << '\n';
			std::cout << std::endl;

			this->recompute_local_segr(active_agent, (move.from == -1) ? move.to : move.from);
			timeline_avg_segr_file.open(this->timeline_file_name + "_avg_segr_timeline.txt", std::ofstream::out | std::ofstream::app);
			timeline_moves_file.open(this->timeline_file_name + "_moves_timeline.txt", std::ofstream::out | std::ofstream::app);

			if (timeline_avg_segr_file.is_open()) {
				timeline_avg_segr_file << std::endl << compute_avg_segr();
				timeline_avg_segr_file.close();
			}
			else {
				std::cout << "cannot open timeline file for output";
				exit(0);
			}

			if (timeline_moves_file.is_open()) {
				if (move.step_name == 'a')
					current_num_of_edges++;
				else if (move.step_name == 'd')
					current_num_of_edges--;
				timeline_moves_file << std::endl << current_num_of_edges << '	' << active_agent << " does " << move << " " << this->graph.color(active_agent) << this->graph.color((move.from == -1) ? move.to : move.from);
				timeline_moves_file.close();
			}
			else {
				std::cout << "cannot open timeline file for output";
				exit(0);
			}
		}
	}
	double cost = social_cost();

	std::string extra_notes = "";
	if (radius_of_steps < graph.graph_num_of_nodes()) {
		extra_notes += "_" + std::to_string(radius_of_steps) + "_local";
	}
	else
		extra_notes += "_global";

	extra_notes += "_" + std::to_string(graph.avg_clustering()) + "_avg_clust_" + std::to_string(this->graph.diameter()) + "_D_" + std::to_string(count) + "_steps";
	this->graph.output_graph_to_dot_file("GE_" + extra_notes_in_output_filename + extra_notes + ".dot");


	//Compute average segregation	
	avg_segr = this->compute_avg_segr();
	std::ofstream out_segr_file(extra_notes_in_output_filename + "_" + std::to_string(avg_segr) + "_segregation.dot");

	if (out_segr_file.is_open())
	{
		for (int i = 0; i < this->graph.graph_num_of_nodes(); i++)
			out_segr_file << local_segregation[i] << '\n';
		out_segr_file.close();
	}
	else
	{
		std::cout << "Cannot open out_file for the segregation output";
		exit(0);
	}

	return avg_segr;
}