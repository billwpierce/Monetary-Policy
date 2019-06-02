#include <math.h>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

struct State {
    double t;
    double u;
    double i;
    double du;
    double di;
    double m;
    double g;
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

int main() {
    double dt = 0.001;
    int num_cycles = (int)ceil(((40 * M_PI) / 23) / dt);
    vector<double> phi_choices = {1, 2, M_PI};
    State state_initial;
    state_initial.t = 0;
    state_initial.u = .039;
    state_initial.i = .024;

    for (int n = 0; n < phi_choices.size(); n++) {
        double phi = phi_choices[n];
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
            state_map.push_back(curr_state);
        }
        double recession_percentage = recession_count / num_cycles;
        double high_unemployment_percentage =
            high_unemployment_count / num_cycles;
        cout << "Phi: " << phi << "\n";
        cout << "Recession Percentage: " << recession_percentage * 100 << "%\n";
        cout << "High Unemployment Percentage: "
             << high_unemployment_percentage * 100 << "%\n";
    }
}