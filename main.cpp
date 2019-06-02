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
};

double calc_du(double i, double t, double phi) {
    return 0.04 * i + 0.012 + 0.0092 * (cos((1.15 * t)) + phi);
}

double calc_di(double u) { return (-0.73 * u) + 0.0438; }

State calc_state(State prev_state, double dt, double phi) {
    State curr_state;
    curr_state.t = prev_state.t + dt;
    curr_state.i = prev_state.i + (dt * prev_state.di);
    curr_state.u = prev_state.u + (dt * prev_state.du);
    curr_state.di = calc_di(curr_state.u);
    curr_state.du = calc_du(curr_state.i, curr_state.t, phi);
    return curr_state;
}

int main() {
    double dt = 0.001;
    vector<double> phi_choices = {1, 2, M_PI};
    State state_initial;
    state_initial.t = 0;
    state_initial.u = .039;
    state_initial.i = .024;

    for (int n = 0; n < phi_choices.size(); n++) {
        double phi = phi_choices[n];
        vector<State> state_map = {state_initial};
        state_map[0].du = calc_du(state_map[0].i, state_map[0].t, phi);
        for (int c = 1; c < 10; c++) {
            State prev_state = state_map[c - 1];
            state_map.push_back(calc_state(prev_state, dt, phi));
        }
    }
}