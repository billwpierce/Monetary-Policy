#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

ofstream outputFile;
std::string filename = "output.csv";

struct State {
    double t;
    double u;
    double i;
    double du;
    double di;
    double m;
    double g;
};

struct State_Results {
    double phi;
    double recession_percentage;
    double high_unemployment_percentage;
};

double calc_du(double i, double t, double phi) {
    return 0.4 * i + 0.012 + 0.0092 * cos(1.15 * t + phi);
}

double calc_di(double u) { return (-0.73 * u) + 0.0438; }

double calc_m(double t, double phi) { return 0.02 * sin(1.15 * t + phi); }

State calc_state(State prev_state, double dt, double phi) {
    State curr_state;
    curr_state.t = prev_state.t + dt;
    curr_state.i = prev_state.i + (dt * prev_state.di);
    curr_state.u = prev_state.u + (dt * prev_state.du);
    curr_state.di = calc_di(curr_state.u);
    curr_state.du = calc_du(curr_state.i, curr_state.t, phi);
    curr_state.m = calc_m(curr_state.t, phi);
    curr_state.g = curr_state.m - curr_state.i;
    return curr_state;
}

State_Results calc_phi_results(State state_initial, double phi,
                               double num_cycles, double dt,
                               bool write_to_csv) {
    if (write_to_csv) {
        outputFile.open(filename);
        outputFile << "U (Unemployment)"
                   << ","
                   << "G (GDP)"
                   << ","
                   << "I (Inflation)"
                   << ","
                   << "M (Money Supply)" << std::endl;
    }
    vector<State> state_map = {state_initial};
    state_map[0].du = calc_du(state_map[0].i, state_map[0].t, phi);
    state_map[0].m = calc_m(state_map[0].t, phi);
    state_map[0].g = state_map[0].m - state_map[0].i;
    double recession_count = 0;
    double high_unemployment_count = 0;
    if (state_map[0].g < 0) {
        recession_count += 1;
    }
    if (state_map[0].u >= .08) {
        high_unemployment_count += 1;
    }
    for (int c = 1; c < num_cycles; c++) {
        State prev_state = state_map[c - 1];
        State curr_state = calc_state(prev_state, dt, phi);
        if (curr_state.g < 0) {
            recession_count += 1;
        }
        if (curr_state.u >= .08) {
            high_unemployment_count += 1;
        }
        if (write_to_csv) {
            outputFile << curr_state.u << "," << curr_state.g << ","
                       << curr_state.i << "," << curr_state.m << std::endl;
        }
        state_map.push_back(curr_state);
    }
    State_Results results;
    results.phi = phi;
    results.recession_percentage = recession_count / num_cycles;
    results.high_unemployment_percentage = high_unemployment_count / num_cycles;
    if (write_to_csv) {
        outputFile.close();
    }
    return results;
}

int main(int argc, char** argv) {
    double dt = 0.001;
    int num_cycles = (int)2 * ceil(((40 * M_PI) / 23) / dt);
    double phi_init = 3;
    double accuracy_digit = 5;
    State state_initial;
    state_initial.t = 0;
    state_initial.u = .039;
    state_initial.i = .024;

    double curr_phi = phi_init;
    double curr_accuracy_digit = 0;
    double run = true;
    while (run) {  // Currently optimizes for unemployment percentage
        State_Results same_phi =
            calc_phi_results(state_initial, curr_phi, num_cycles, dt, false);
        cout << "Phi: " << curr_phi << "\n";
        cout << "Recession Percentage: " << same_phi.recession_percentage * 100
             << "%\n";
        cout << "High Unemployment Percentage: "
             << same_phi.high_unemployment_percentage * 100 << "%\n";
        double high_phi = curr_phi + pow(10, -1 * curr_accuracy_digit);
        double low_phi = curr_phi - pow(10, -1 * curr_accuracy_digit);
        State_Results hp_results =
            calc_phi_results(state_initial, high_phi, num_cycles, dt, false);
        State_Results lp_results =
            calc_phi_results(state_initial, low_phi, num_cycles, dt, false);
        if (hp_results.recession_percentage < same_phi.recession_percentage) {
            if (lp_results.recession_percentage <
                hp_results.recession_percentage) {
                curr_phi = low_phi;
            } else {
                curr_phi = high_phi;
            }
        } else if (lp_results.recession_percentage <
                   same_phi.recession_percentage) {
            curr_phi = low_phi;
        } else {
            if (curr_accuracy_digit >= accuracy_digit) {
                run = false;
            } else {
                curr_accuracy_digit += 1;
            }
        }
    }
    calc_phi_results(state_initial, curr_phi, num_cycles, dt, false);
    return 0;
}