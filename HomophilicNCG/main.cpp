#include "graph.h"
#include "NCG.h"

//#include<math.h>

//for Windows
//#include<conio.h>

int my_log(int base, int value) {
	double result = log(value) / log(base);
	return static_cast<int>(result);
}

int main() {
	srand(time(0));

//extra parameters. Don't change for this model
	int infty = 100000;
	bool round_robin = 0;
	bool storage_distances = 1;
	bool rand_endpoint = 1;

	//global variables. do not change
	std::string extra_notes;
	std::string if_pairwise;// = (pairwise ? "PSN" : "NE");
	std::string if_BR;// = (best_resp ? "BR" : "IR");

	//to define
	bool if_rand_color= 0; //random color distribution; otherwise, color is assigned in a greedy way, i.e., one color to first num_of_1st_color_nodes nodes
	bool best_resp = 1;
	std::string model_name = "ICFNCG"; //DEINCG or ICFNCG
	double apx_factor = 1.0; //to generate approximately stable networks
	std::string init_graph = "tree";//"tree" or "grid"
	bool pairwise = 1; 
	std::set<char> list_of_prohib_moves = {  'd' }; // 'd' means no deletions (add-only version), 'a' for no edge additions

		for (auto& n : { 1000 }) {
			for (auto& num_of_1st_color_nodes : { 500 })
			{
				for (double alpha = 10.0; alpha <= 25.0; alpha += 5.0)
				{
					//create full name of the configuration
					if_BR = (best_resp ? "BR" : "IR");
					std::string graph_name = (if_rand_color) ? ("rand_" + init_graph) : ("fixed_" + init_graph);

					std::string concept =(pairwise)? "PSN":"";
					
					if (list_of_prohib_moves.find('d') == list_of_prohib_moves.end()) 
					{
						if (apx_factor > 1.0) {
							std::string apx_val = std::to_string(apx_factor);
							apx_val.erase(apx_val.find_last_not_of('0') + 1, std::string::npos);
							concept = apx_val + concept;
						}
					}
					else
					{
						concept = "Ao" + concept;
						std::string apx_val = std::to_string(apx_factor);
						apx_val.erase(apx_val.find_last_not_of('0') + 1, std::string::npos);
						if (apx_factor > 1.0)concept = apx_val + concept;
					}

					std::string alpha_val = std::to_string((int)alpha);

					extra_notes = std::to_string(n) + "_n_" + alpha_val + "_alpha_" + std::to_string(num_of_1st_color_nodes) + "_of_black_" 
								+ model_name+"_"+concept+"_from_" + graph_name+"_rand_" + if_BR;
					//***************


					std::ofstream avg_segr_file(extra_notes + "_avg_segregation.dot");

					//one round correspond to one run of the game, one run ends when the game converges to a stable state or approximately stable state 
					for (int num_of_round = 1; num_of_round <= 50; num_of_round++)
					{
						Node_colored_graph g_color(n, infty, num_of_1st_color_nodes, if_rand_color);

						if(init_graph == "tree")
							g_color.tree_graph(2);
						else if (init_graph=="grid")
							g_color.grid_graph(50, 20);
						else {
							std::cout << "initial grap is not defined";
							exit(0);
						}

						std::string file_name_for_one_round;

						file_name_for_one_round = extra_notes + "_" + std::to_string(num_of_round) + "_run";

						SchellingNCG game_apx(g_color, alpha, list_of_prohib_moves, file_name_for_one_round);
						double avg_segr = game_apx.greedy_find_GE(file_name_for_one_round, round_robin, rand_endpoint, best_resp, -1, 0, apx_factor);

						if (avg_segr_file.is_open()) {
							avg_segr_file << avg_segr << std::endl;
						}
						else {
							std::cout << "Cannot open out_file for the segregation output";
							exit(0);
						}
					}
				}
			}

		}

	std::cout << "\n done \n";
	return 0;
}
 