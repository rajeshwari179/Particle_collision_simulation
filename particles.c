#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Define a structure to represent a particle
typedef struct
{
    double x, y, x_velocity, y_velocity, lastup_time;
    int w_collision, p_collision, label;
} particle;
// Define a structure to represent a collision
typedef struct
{
    particle particle1, particle2;
    double collision_time;
} collision;
// Function to calculate the time until a wall collision occurs for a particle
particle wall_collision_time(particle *p, double Lx, double Ly, double R)
{
    // Calculate potential collision times with each wall
    double leftedge_col = p->lastup_time + (R - p->x) / p->x_velocity;
    double rightedge_col = p->lastup_time + ((Lx - R) - p->x) / p->x_velocity;
    double bottom_col = p->lastup_time + (R - p->y) / p->y_velocity;
    double top_col = p->lastup_time + ((Ly - R) - p->y) / p->y_velocity;
    if (p->x_velocity > 0)
    {
        leftedge_col = 1e10; // If the particle is moving away from a wall, set the corresponding collision time to a large value
    }
    else
    {
        rightedge_col = 1e10;
    }
    if (p->y_velocity > 0)
    {
        bottom_col = 1e10;
    }
    else
    {
        top_col = 1e10;
    }
    double min_value = 1e10; // Initialize to a large positive value
    int ref_value = 0;       // Variable to store the reference value

    if (leftedge_col >= 0 && leftedge_col < min_value)
    {
        min_value = leftedge_col;
        ref_value = -1; // left wall label
    }

    if (rightedge_col >= 0 && rightedge_col < min_value)
    {
        min_value = rightedge_col;
        ref_value = -2; // right wall label
    }

    if (bottom_col >= 0 && bottom_col < min_value)
    {
        min_value = bottom_col;
        ref_value = -3; // bottom wall label
    }

    if (top_col >= 0 && top_col < min_value)
    {
        min_value = top_col;
        ref_value = -4; // top wall label
    }

    particle p1;
    p1.label = ref_value;
    p1.lastup_time = min_value;
    return p1; // funtion returns the particle which contains time required for collision in the variable lastup_time and wall label
}
// Function to calculate the time until a particle-particle collision occurs

double particle_collision_time(particle *p1, particle *p2, double R)
{

    double time = fmax(p1->lastup_time, p2->lastup_time); // take the maximum of lasst updated time of two particles
    double dummyp2_x, dummyp2_y, dummyp1_x, dummyp1_y;
    // bring the both particles to the same time
    if (p1->lastup_time > p2->lastup_time)
    {
        dummyp2_x = p2->x + (p2->x_velocity) * (p1->lastup_time - p2->lastup_time);
        dummyp2_y = p2->y + (p2->y_velocity) * (p1->lastup_time - p2->lastup_time);
        dummyp1_x = p1->x;
        dummyp1_y = p1->y;
    }
    else
    {
        dummyp1_x = p1->x + (p1->x_velocity) * (p2->lastup_time - p1->lastup_time);
        dummyp1_y = p1->y + (p1->y_velocity) * (p2->lastup_time - p1->lastup_time);
        dummyp2_x = p2->x;
        dummyp2_y = p2->y;
    }
    // solving the quadratic equation
    double dx = dummyp1_x - dummyp2_x;
    double dy = dummyp1_y - dummyp2_y;
    double dvx = p1->x_velocity - p2->x_velocity;
    double dvy = p1->y_velocity - p2->y_velocity;
    double a = dvx * dvx + dvy * dvy;
    double b = 2 * (dx * dvx + dy * dvy);
    double c = dx * dx + dy * dy - 4 * R * R;

    double D = b * b - 4 * a * c;
    // discriminant negative set the collision time to infinity

    if (D < 0)
    {
        return 1e10;
    }
    // discriminant = 0, the equation has one root

    else if (D == 0)
    {
        double t = (-b + sqrt(D)) / (2 * a);
        if (t > 0)
        {
            return t + time;
        }
    }
    else // discriminant greater than 0, two roots  take the minimum positive root

    {
        double t1 = (-b - sqrt(D)) / (2 * a);
        double t2 = (-b + sqrt(D)) / (2 * a);

        if (t1 > 0 && t2 > 0)
        {
            return fmin(t1, t2) + time;
        }
        else if (t1 > 0)
        {
            return t1 + time;
        }
        else if (t2 > 0)
        {
            return t2 + time;
        }
    }
    return 1e10;
}
// Function to perform insertion sort on an array of collisions

void insertion_sort(collision array[], int N)
{
    int i, j;
    collision temp;
    int n = N * (N + 1) / 2;
    for (i = 1; i < n; i++)
    {
        temp = array[i];
        j = i - 1;
        while (j >= 0 && temp.collision_time < array[j].collision_time)
        {
            array[j + 1] = array[j];
            j -= 1;
        }
        array[j + 1] = temp;
    }
}
// Function to update particle positions and velocities after a collision

void update(particle particles[], collision *collision, int N)
{
    // store the labels of particles involved in collision

    int idx1 = collision->particle1.label;
    int idx2 = collision->particle2.label;
    // update the position of the particles involved in collision

    particles[idx1].x += (particles[idx1].x_velocity) * (collision->collision_time - collision->particle1.lastup_time);

    particles[idx1].y += (particles[idx1].y_velocity) * (collision->collision_time - collision->particle1.lastup_time);

    particles[idx1].lastup_time = collision->collision_time; // update the last updated time of particles

    if (idx2 >= 0) // check if the collision is particle-particle
    {
        particles[idx2].x += (particles[idx2].x_velocity) * (collision->collision_time - collision->particle2.lastup_time);
        particles[idx2].y += (particles[idx2].y_velocity) * (collision->collision_time - collision->particle2.lastup_time);
        particles[idx2].lastup_time = collision->collision_time;
        particles[idx1].p_collision++;
        particles[idx2].p_collision++;
    }

    if (idx2 < 0) // check if the collision is particle-wall
    {
        particles[idx1].w_collision++;
        particles[idx2].w_collision++;

        if (idx2 > -3) // check if the wall is either left or right
        {
            particles[idx1].x_velocity = -particles[idx1].x_velocity;
        }
        else // case where the wall is top or bottom
        {
            particles[idx1].y_velocity = -particles[idx1].y_velocity;
        };
    }
    else // velocity changes if the collision is particle-particle

    {
        // interchange the velocity of particles involved in the collision

        double u1x = particles[idx1].x_velocity;
        double u1y = particles[idx1].y_velocity;
        double u2x = particles[idx2].x_velocity;
        double u2y = particles[idx2].y_velocity;

        particles[idx1].x_velocity = u2x;
        particles[idx1].y_velocity = u2y;

        particles[idx2].x_velocity = u1x;
        particles[idx2].y_velocity = u1y;
    }
}

int main(int argc, char *argv[])
{

    double curr_time = 0; // initialise current time
    char *input_file = argv[1];
    double end_time = atof(argv[2]); // user defined end time

    FILE *file = fopen(input_file, "r");

    int N;            // no of particles
    double R, Lx, Ly; // radius and dimenssions of the box

    // Read values from the input file
    fscanf(file, "%d", &N);
    // printf("%d", N);
    fscanf(file, "%lf", &R);
    // printf("%lf", R);

    fscanf(file, "%lf %lf", &Lx, &Ly);

    // Allocate memory for particles and potential collisions
    particle *particles = malloc(N * sizeof(particle));
    
    collision *collisions = malloc((N * (N + 1) / 2) * sizeof(collision));

    // Check if memory allocation was successful
    if (particles == NULL)
    {
        perror("Memory allocation failed");
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        fscanf(file, "%lf %lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].x_velocity, &particles[i].y_velocity);

        particles[i].label = i;
        particles[i].lastup_time = 0.0;
        particles[i].w_collision = 0;
        particles[i].p_collision = 0;
    }
    int l1 = -5;
    int l2 = -5;

    while (curr_time < end_time) // stopping condition of the process

    {

        int z = 0;
        // append particle-particle collision to the array of collisions

        for (int i = 0; i < N - 1; i++)
        {

            for (int j = i + 1; j < N; j++)
            {

                double temptime = particle_collision_time(&particles[i], &particles[j], R);
                // condition to check whether particles in the collsions are not same as the previous collision

                if (particles[i].label != l1 && particles[j].label != l2 && temptime != curr_time || curr_time == 0)
                {
                    collisions[z].collision_time = particle_collision_time(&particles[i], &particles[j], R);

                    collisions[z].particle1 = particles[i];
                    collisions[z].particle2 = particles[j];
                    z++;
                }
            }
        }
        // append particle-wall collision to the array of collisions

        for (int c = 0; c < N; c++)
        {

            particle p3 = wall_collision_time(&particles[c], Lx, Ly, R);

            collisions[z].particle1 = particles[c];
            collisions[z].collision_time = p3.lastup_time;
            collisions[z].particle2 = p3;

            z++;
        }

        insertion_sort(collisions, N); // sort the order of collisions wrt time

        if (collisions[0].collision_time <= end_time)
        {
            update(particles, &collisions[0], N);
        }
        // double temp_time = curr_time;

        int p = 1;
        double temp_time = curr_time;
        curr_time = collisions[0].collision_time;
        l1 = collisions[0].particle1.label;
        l2 = collisions[0].particle2.label;

        curr_time = collisions[0].collision_time;
    }
    
    // bring all particles to the endtime position
    for (int i = 0; i < N; i++)
    {

        particles[i].x += (particles[i].x_velocity) * (-particles[i].lastup_time + end_time);
        particles[i].y += (particles[i].y_velocity) * (-particles[i].lastup_time + end_time);
        printf("%.6f, %.6f, %d, %d\n", particles[i].x, particles[i].y, particles[i].w_collision, particles[i].p_collision);
        
    }

    fclose(file);
    free(collisions);
    free(particles);
}