#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include "ortools/sat/cp_model.h"

using namespace std;
namespace operations_research::sat{
    struct Instance{
        int num_locations;
        vector<vector<int>> distance;
        ifstream fin;
        string type_data;
        string file_name;

        Instance(string file_name):file_name(file_name){
            fin = ifstream(file_name);
            fin >> num_locations;
            fin >> type_data;
            distance = vector<vector<int>>(num_locations, vector<int>(num_locations, 0));
            if(type_data == "EUC_2D"){
                read_distance_data_euc2D_coord_matrix();
            }else if(type_data == "ATT"){
                read_distance_data_att_coord_matrix();
            }else if(type_data == "GEO"){
                read_distance_data_geo_coord_matrix();
            }else if(type_data == "LOWER_DIAG_ROW"){
                read_distance_data_lower_diag_matrix();
            }else if(type_data == "UPPER_DIAG_ROW"){
                read_distance_data_upper_diag_matrix();
            }else if(type_data == "FULL_MATRIX"){
                read_distance_data_full_matrix();
            }else if(type_data == "CEIL_2D"){
                read_distance_data_ceil2D_coord_matrix();
            }else if(type_data == "UPPER_ROW"){
                read_distance_data_upper_matrix();
            }else if(type_data == "LOWER_ROW"){
                read_distance_data_lower_matrix();
            }else{
                cout << "Type of data is not supported\n";
            }
            cout << "-------------------Read data successfully!-------------------\n";
        }

        void read_distance_data_full_matrix(){
            for (int i = 0; i < num_locations; i++)
                for (int j = 0; j < num_locations; j++) {
                    fin >> distance[i][j];
                }
        }

        void read_distance_data_upper_matrix(){
            for (int i = 0; i < num_locations; i++) {
                for (int j = i + 1; j < num_locations; j++) {
                    fin >> distance[i][j];
                }
                for (int j = 0; j < i; j++) {
                    distance[i][j] = distance[j][i];
                }
            }
        }

        void read_distance_data_lower_matrix(){
            double tmp_number;
            for (int i = 1; i < num_locations; i++)
                for (int j = 0; j < i; j++) {
                    fin >> tmp_number;
                    distance[i][j] = distance[j][i] = tmp_number;
                }
        }

        void read_distance_data_geo_coord_matrix(){
            double temp;
            vector<double> x;
            vector<double> y;
            vector<double> latitude;
            vector<double> longitude;
            for (int i = 0; i < num_locations; i++) {
                x.push_back(0);
                y.push_back(0);
                latitude.push_back(0);
                longitude.push_back(0);
            }
            for (int i = 0; i < num_locations; i++) {
                fin >> temp;
                fin >> x[i];
                fin >> y[i];
            }
            double PI = 3.141592, RRR = 6378.388;
            for (int i = 0; i < num_locations; i++){
                double deg, min;
                //deg = (int) (x[i]+0.5);
                deg = (int) x[i];
                min = x[i] - deg;
                latitude[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
                //deg = (int) (y[i]+0.5);
                deg = (int) y[i];
                min = y[i] - deg;
                longitude[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
            }
            for (int i = 0; i < num_locations; i++)
                for (int j = 0; j < num_locations; j++) {
                    double q1 = cos(longitude[i] - longitude[j]);
                    double q2 = cos(latitude[i] - latitude[j]);
                    double q3 = cos(latitude[i] + latitude[j]);
                    distance[i][j] = (i != j) ? (int) (RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0) : 0;
                }

        }

        void read_distance_data_att_coord_matrix(){
            double temp;
            vector<double> x;
            vector<double> y;
            for (int i = 0; i < num_locations; i++) {
                x.push_back(0);
                y.push_back(0);
            }
            for (int i = 0; i < num_locations; i++) {
                fin >> temp;
                fin >> x[i];
                fin >> y[i];
            }
            for (int i  = 0; i < num_locations; i++)
                for (int j = 0; j < num_locations; j++){
                    double xd = x[i] - x[j];
                    double yd = y[i] - y[j];
                    double rij = sqrt((xd*xd + yd*yd) / 10.0);
                    double tij = (int) (rij+0.5);
                    distance[i][j] = tij < rij ? tij + 1 : tij;
                }
        }

        void read_distance_data_ceil2D_coord_matrix(){
            double temp;
            vector<double> x;
            vector<double> y;
            for (int i = 0; i < num_locations; i++) {
                x.push_back(0);
                y.push_back(0);
            }
            for (int i = 0; i < num_locations; i++) {
                fin >> temp;
                fin >> x[i];
                fin >> y[i];
            }
            for (int i = 0; i < num_locations; i++)
                for (int j = 0; j < num_locations; j++) {
                    double xd = x[i] - x[j];
                    double yd = y[i] - y[j];
                    distance[i][j] = (int)(sqrt(xd * xd + yd * yd) + 0.5);
                }
        }

        void read_distance_data_euc2D_coord_matrix(){
            double temp;
            vector<double> x;
            vector<double> y;
            for (int i = 0; i < num_locations; i++) {
                x.push_back(0);
                y.push_back(0);
            }
            for (int i = 0; i < num_locations; i++) {
                fin >> temp;
                fin >> x[i];
                fin >> y[i];
            }
            for (int i  = 0; i < num_locations; i++)
                for (int j = 0; j < num_locations; j++){
                    double xd = x[i] - x[j];
                    double yd = y[i] - y[j];
                    distance[i][j] = (int)(sqrt(xd*xd + yd*yd)+0.5);
                }
        }

        void read_distance_data_lower_diag_matrix() {
            double tmp_number;
            for (int i = 0; i < num_locations; i++) {
                for (int j = 0; j <= i; j++) {
                    fin >> tmp_number;
                    distance[i][j] = distance[j][i] = tmp_number;
                }
            }
        }

        void read_distance_data_upper_diag_matrix() {
            for (int i = 0; i < num_locations; i++) {
                for (int j = i; j < num_locations; j++) {
                    fin >> distance[i][j];
                    distance[j][i] = distance[i][j];
                }
            }
        }
    };

    void printData(const Instance &instance){
        cout << "Number of locations: " << instance.num_locations << endl;
        cout << "Distance:" << endl;
        cout << "    ";
        for (int i = 0; i < instance.num_locations; i++)
            printf("%4d ", i);
        cout << endl;
        cout << "   ";
        for (int i = 0; i < instance.num_locations; i++)
            printf("_____");
        cout << endl;
        for (int i = 0; i < instance.num_locations; i++) {
            printf("%2d| ", i);
            for (int j = 0; j < instance.num_locations; j++) {
                printf("%4d ", instance.distance[i][j]);
            }
            cout << endl;
        }
    }

    void SolveTSP_AllDifferent(const Instance &instance, string outfile){
        vector<vector<int64_t>> d(instance.num_locations, vector<int64_t>(instance.num_locations));
        for(int i = 0; i < instance.num_locations; ++i)
            for(int j = 0; j < instance.num_locations; ++j)
                d[i][j] = (int64_t) instance.distance[i][j];

        CpModelBuilder cp_model;

        int max = 0;
        for(int i = 0; i < instance.num_locations; ++i)
            for(int j = 0; j < instance.num_locations; ++j)
                if (max < d[i][j])
                    max = d[i][j];

        const Domain domain(0, instance.num_locations-1);
        const Domain domainCost(1, max);

        vector<IntVar> x, next, prev, z;
        x.reserve(instance.num_locations);
        next.reserve(instance.num_locations);
        prev.reserve(instance.num_locations);
        z.reserve(instance.num_locations);
        for(int i = 0; i < instance.num_locations; ++i) {
            x.push_back(cp_model.NewIntVar(domain));
            next.push_back(cp_model.NewIntVar(domain));
            prev.push_back(cp_model.NewIntVar(domain));
            z.push_back(cp_model.NewIntVar(domainCost));
        }
        cp_model.AddEquality(x[0], 0);
        cp_model.AddAllDifferent(x);
        cp_model.AddAllDifferent(next);
        for(int i = 1; i < instance.num_locations; ++i)
            cp_model.AddVariableElement(x[i-1], next, x[i]);
        cp_model.AddAllDifferent(prev);
        cp_model.AddInverseConstraint(prev, next);
        cp_model.AddLessThan(next[0], prev[0]);

        for (int i = 0; i < instance.num_locations; ++i)
            cp_model.AddElement(next[i], d[i], z[i]);

        LinearExpr obj;
        for (int i = 0; i < instance.num_locations; ++i)
            obj += z[i];
        cp_model.Minimize(obj);

        SatParameters parameters;
        parameters.set_log_search_progress(true);
        parameters.set_max_time_in_seconds(3600);
        parameters.set_use_combined_no_overlap(true);
        parameters.set_num_search_workers(1);

        Model model;
        model.Add(NewSatParameters(parameters));

        const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);

        if (response.status() == CpSolverStatus::OPTIMAL ||
            response.status() == CpSolverStatus::FEASIBLE) {
            cout << "Best bound: " << response.best_objective_bound() << endl;
            cout << "Lower bound: " << response.inner_objective_lower_bound() << endl;
            cout << "Gap: " << response.gap_integral() << endl;
            cout << "Objective: " << response.objective_value() << endl;
            cout << "Time: " << response.user_time() << endl;
            vector<int> solx(instance.num_locations);
            vector<int> solz(instance.num_locations);
            for(int i = 0; i < instance.num_locations; ++i){
                solx[i] = SolutionIntegerValue(response, x[i]);
                solz[i] = SolutionIntegerValue(response, z[i]);
            }
            cout << "Solution: \n";
            for (int i = 0; i < solx.size(); i++)
                printf("%d ", solx[i]);
            cout << "0";
            int distance = 0;
            for (int i = 0; i < solz.size(); i++){
                distance += solz[i];
            }
            cout << "\nDistance: " << distance << endl;
            //LOG(INFO) << response.best_objective_bound() << " " << response.objective_value() << " " << response.user_time();
            ofstream fout;
            fout.open(outfile);
            fout << "Objective: " << response.objective_value() << endl;
            fout << "Lower bound: " << response.best_objective_bound() << endl;
            fout << "Time: " << response.user_time() << endl;
            fout << "Solution: ";
            for (int i = 0; i < solx.size(); i++)
                fout << solx[i] << " ";
            fout << "0";
        } else {
            LOG(INFO) << "No solution found.";
        }
    }

    void SolveTSP_Circuit(const Instance &instance, string outfile){
        CpModelBuilder cp_model;
        vector<vector<BoolVar>> x(instance.num_locations, vector<BoolVar>(instance.num_locations));
        LinearExpr obj;
        auto circuit = cp_model.AddCircuitConstraint();
        for(int i = 0; i < instance.num_locations; ++i)
            for(int j = 0; j < instance.num_locations; ++j){
                x[i][j] = cp_model.NewBoolVar();
                if(i != j) {
                    obj += x[i][j] * instance.distance[i][j];
                    circuit.AddArc(i, j, x[i][j]);
                }
            }
        cp_model.Minimize(obj);

        SatParameters parameters;
        parameters.set_log_search_progress(true);
        parameters.set_max_time_in_seconds(3600);
        parameters.set_use_combined_no_overlap(true);
        parameters.set_num_search_workers(1);
        Model model;

        model.Add(NewSatParameters(parameters));

        const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);

        if (response.status() == CpSolverStatus::OPTIMAL ||
            response.status() == CpSolverStatus::FEASIBLE) {
            cout << "Best bound: " << response.best_objective_bound() << endl;
            cout << "Lower bound: " << response.inner_objective_lower_bound() << endl;
            cout << "Gap: " << response.gap_integral() << endl;
            cout << "Objective: " << response.objective_value() << endl;
            cout << "Time: " << response.user_time() << endl;
            vector<vector<int>> sol_x(instance.num_locations, vector<int>(instance.num_locations));
            vector<int> tour;
            for(int i = 0; i < instance.num_locations; ++i)
                for(int j = 0; j < instance.num_locations; ++j){
                    if(i != j) {
                        sol_x[i][j] = SolutionIntegerValue(response, x[i][j]);
                    }
                }

            tour.push_back(0);
            while(tour.size() < instance.num_locations){
                for(int i = 1; i < instance.num_locations; i++)
                    if(sol_x[tour[tour.size()-1]][i] == 1){
                        tour.push_back(i);
                    }
            }
            tour.push_back(0);
            cout << "Solution: ";
            for(int i = 0; i < tour.size(); i++)
                cout << tour[i] << " ";
            //LOG(INFO) << response.best_objective_bound() << " " << response.objective_value() << " " << response.user_time();
            ofstream fout;
            fout.open(outfile);
            fout << "Objective: " << response.objective_value() << endl;
            fout << "Lower bound: " << response.best_objective_bound() << endl;
            fout << "Time: " << response.user_time() << endl;
            fout << "Solution: ";
            for(int i = 0; i < tour.size(); i++)
                fout << tour[i] << " ";
        } else {
            LOG(INFO) << "No solution found.";
        }
    }
}

vector<string> SplitStringWithDelimiter(string s, string delimiter) {
    vector<string> returnValue;
    string::size_type start = 0;
    string::size_type end = s.find(delimiter);

    while(end != string::npos) {
        returnValue.push_back(s.substr(start, end-start));
        start = end + 1;
        end = s.find(delimiter, start);
    }

    returnValue.push_back(s.substr(start));
    return returnValue;
}


int main(int argc, char* argv[]){
    string data = argv[1];
    string datafile = data + ".tsp";
    cout << "Instance name: " << datafile << endl;
    auto fn = SplitStringWithDelimiter(data, "/");
    string outfile = "/home/cuong/CLionProjects/TSP_CP/tsp_ortools_result/" + fn[fn.size()-1] + ".sol";
    if (filesystem::exists(outfile)) {
        cout << "Already solved!!!: " << datafile << endl;
        return EXIT_SUCCESS;
    }
//    string model = argv[2];
    operations_research::sat::Instance instance(data);
    operations_research::sat::printData(instance);
    operations_research::sat::SolveTSP_Circuit(instance, outfile);
//    if(model == "diff"){
//        string outfile = data + "_diff.sol";
//        operations_research::sat::SolveTSP_AllDifferent(instance, outfile);
//    }
//    if(model == "circuit") {
//        string outfile = data + "_circuit.sol";
//        operations_research::sat::SolveTSP_Circuit(instance, outfile);
//    }
    return EXIT_SUCCESS;
}