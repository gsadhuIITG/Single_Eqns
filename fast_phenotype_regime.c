#include <stdio.h>
#include <math.h>

#define NX 1001        // spatial grid points
#define NT 5000        // time steps

int main() {

    int i, n;

    /* ================= Parameters ================= */
    double L = 10.0;
    double T = 50.0;
    double K = 1.0;

    /* Population & switching */
    double mu_p      = 1.0;
    double lambda_pm = 20.0;   // fast switching
    double lambda_mp = 20.0;
    double Dm        = 0.001;

    /* Hill function */
    double cH = 0.5;
    int    k  = 3;

    /* Cyclic hypoxia */
    double c1 = 1.0;      // normoxia (> cH)
    double c2 = 0.2;      // hypoxia   (< cH)
    double Tcycle = 10.0;

    /* ================= Discretization ================= */
    double dx = L / (NX - 1);
    double dt = T / NT;

    /* ================= Arrays ================= */
    double p[NX], m[NX], s[NX];
    double p_new[NX], m_new[NX], s_new[NX];

    /* ================= Initial conditions ================= */
    for (i = 0; i < NX; i++) {
        double x = i * dx;

        if (x <= 2.0) {
            p[i] = 0.2;
            s[i] = 0.2;
        } else {
            p[i] = 0.0;
            s[i] = 0.0;
        }

        m[i] = 0.0;
    }

    /* ================= Output control ================= */
    int output_index = 0;
    double next_output_time = 0.0;

    /* ================= Time loop ================= */
    for (n = 0; n <= NT; n++) {

        double t = n * dt;

        /* ===== Cyclic hypoxia ===== */
        double c_t;
        if (fmod(t, Tcycle) < 0.5 * Tcycle)
            c_t = c1;
        else
            c_t = c2;

        /* Hill + rates */
        double Psi = pow(c_t, k) / (pow(cH, k) + pow(c_t, k));
        double mu      = mu_p * Psi;
        double lambda1 = lambda_pm * (1.0 - Psi);
        double lambda2 = lambda_mp * Psi;
        double Theta_c = lambda1 / (lambda1 + lambda2);

        /* ===== OUTPUT (single correct block) ===== */
        if (fabs(t - next_output_time) < 0.5 * dt && output_index <= 50) {

            char fname_pm[100], fname_s[100];
            sprintf(fname_pm, "pm_t_%d.txt", output_index);
            sprintf(fname_s,  "s_t_%d.txt",  output_index);

            FILE *fp = fopen(fname_pm, "w");
            FILE *fs = fopen(fname_s,  "w");

            for (i = 0; i < NX; i++) {
                fprintf(fp, "%0.5f\t%0.8f\t%0.8f\n",
                        i * dx, p[i], m[i]);
                fprintf(fs, "%0.5f\t%0.8f\n",
                        i * dx, s[i]);
            }

            fclose(fp);
            fclose(fs);

            printf("Saved p, m, s at t = %d\n", output_index);

            output_index++;
            next_output_time += 1.0;
        }

        if (n == NT) break;

        /* ===== p equation ===== */
        for (i = 0; i < NX; i++) {
            double total = p[i] + m[i];
            p_new[i] = p[i] + dt * (
                mu * p[i] * (1.0 - total / K)
                - lambda1 * p[i]
                + lambda2 * m[i]
            );
            if (p_new[i] < 0) p_new[i] = 0;
        }

        /* ===== m equation ===== */
        for (i = 1; i < NX - 1; i++) {
            double lap_m = (m[i+1] - 2.0*m[i] + m[i-1]) / (dx * dx);
            m_new[i] = m[i] + dt * (
                Dm * lap_m
                + lambda1 * p[i]
                - lambda2 * m[i]
            );
            if (m_new[i] < 0) m_new[i] = 0;
        }

        m_new[0]      = m_new[1];
        m_new[NX - 1] = m_new[NX - 2];

        for (i = 0; i < NX; i++) {
            p[i] = p_new[i];
            m[i] = m_new[i];
        }

        /* ===== Reduced equation for s ===== */
        for (i = 1; i < NX - 1; i++) {
            double lap_s = (s[i+1] - 2.0*s[i] + s[i-1]) / (dx * dx);
            s_new[i] = s[i] + dt * (
                Dm * Theta_c * lap_s
                + (1.0 - Theta_c) * mu * s[i] * (1.0 - s[i] / K)
            );
            if (s_new[i] < 0) s_new[i] = 0;
        }

        s_new[0]      = s_new[1];
        s_new[NX - 1] = s_new[NX - 2];

        for (i = 0; i < NX; i++)
            s[i] = s_new[i];
    }

    printf("Simulation finished successfully.\n");
    return 0;
}
