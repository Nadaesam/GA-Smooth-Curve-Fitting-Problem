#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>

using namespace std ;

int upper = 10 , lower = -10 ;
double ratio = 0.2 ;
const int population = 6 ;
const int l = ceil(ratio*population) ;
vector<int> mutationLocation ;
const int MaxGeneration = 10 ;
vector<vector<double>> offsprings ;
vector<pair<double , double>> points ;
pair<vector<double>,double> BestSolution ;
void Mutation( vector<vector<double>> &OS , int d , vector<int> mutationLocation , int curGeneration)
{
    float Pm = 0.1 ;
    for (int i = 0; i < mutationLocation.size(); i++)
    {
        double y ;
        for (int j = 0; j < d; j++)
        {
            float randomNum = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if( randomNum <= Pm )
            {
                double deltaL = OS[mutationLocation[i]][j] - lower ;
                double deltaU = upper - OS[mutationLocation[i]][j]  ;
                float randomNum1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                if( randomNum1 <= 0.5)
                {
                    y = deltaL ;
                }
                else
                {
                    y = deltaU ;
                }
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                double b = ((double) rand() / RAND_MAX) * (5 - 1) + 1 ;
                double equation = y * (1 - pow(r,pow( 1 - curGeneration / MaxGeneration , b )))  ;
                if( y == deltaL )
                {
                    OS[mutationLocation[i]][j] -= equation ;
                }
                else
                {
                    OS[mutationLocation[i]][j] += equation ;
                }

            }
        }

    }
}
void CrossOver( vector<vector<double>> &OS , int d )
{
    int i = 0 ;
    float Pc = 0.7 ;
    for ( i = l; i < population-1; i+=2)
    {
        float randomNum = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if( randomNum <= Pc )
        {
            mutationLocation.push_back(i) ;
            mutationLocation.push_back(i+1) ;
            //srand(time(0)) ;
            int random1 = rand()%(d-1)+1 ;
            int random2 = rand()%(d-1)+1 ;
            while( random1 == random2 )
            {
                random1 = rand()%(d-1)+1 ;
            }
            int Min = min(random1,random2) ;
            int Max = max(random1,random2) ;
            for (int j = Min; j < Max; j++)
            {
                swap(OS[i][j] , OS[i+1][j]) ;
            }


        }
    }
    if( i == population-1 )
    {
        float randomNum = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if( randomNum <= Pc )
        {
            mutationLocation.push_back(i) ;
            srand(time(0)) ;
            int random1 = rand()%(d-1)+1 ;
            int random2 = rand()%(d-1)+1 ;
            while( random1 == random2 )
            {
                random1 = rand()%(d-1)+1 ;
            }
            int crossOver = l + (rand() % (population-l));
            int Min = min(random1,random2) ;
            int Max = max(random1,random2) ;
            for (int j = Min; j < Max; j++)
            {
                swap(OS[i][j] , OS[crossOver][j]) ;
            }
        }
    }
}
bool cmp( pair<vector<double>,double> &p1 , pair<vector<double>,double> &p2 )
{
    return p1.second > p2.second ;
}
void tournamentSelection( vector<pair<vector<double>,double>> &sol_fitness  )
{
    sort(sol_fitness.begin(),sol_fitness.end(),cmp) ;
    for (int i = 0; i < l; i++)
    {
        offsprings.push_back(sol_fitness[i].first) ;
    }
    for (int i = l; i < population; i++)
    {
        int index1 = l + (rand() % (population-l));
        int index2 = l + (rand() % (population-l)) ;
        while( index1 == index2 )
        {
            index1 = l + (rand() % (population-l));
        }
        if( sol_fitness[index1].second >= sol_fitness[index2].second )
        {
            offsprings.push_back(sol_fitness[index1].first) ;
        }
        else
        {
            offsprings.push_back(sol_fitness[index2].first) ;
        }

    }
}
vector<double> GenerateRandomSolutions(int degree )
{
    vector<double> solution ;
    for (int i = 0; i < degree; i++)
    {
        double gene = ((double) rand() / RAND_MAX) * (upper - lower) + lower;
        gene = round(gene * 100.0) / 100.0;
        solution.push_back(gene) ;
    }
    return solution ;

}
double meanSquareError( vector<double> solution , int d )
{
    double result = 0 ;
    for (int i = 0; i < points.size(); i++)
    {
        double subresult = 0 ;
        for (int j = 0; j < d; j++)
        {
            if( points[i].first == 0.0 )
            {
                subresult += solution[j] ;
            }
            else
            {
                subresult += solution[j] * pow(points[i].first,j) ;
            }
        }
        result += pow(subresult-points[i].second , 2) ;

    }
    return result/points.size() ;
}
double fitness( vector<double> solution , int d )
{
    double result = meanSquareError(solution,d) ;

    return 1/result ;
}
int main() {
    srand(time(0));

    ifstream inputFile("curve_fitting_input.txt");
    ofstream outputFile("output.txt");

    if (!inputFile.is_open()) {
        cerr << "Error opening input file." << endl;
        return 1;
    }

    if (!outputFile.is_open()) {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    int numDatasets;
    inputFile >> numDatasets;

    for (int datasetIndex = 0; datasetIndex < numDatasets; ++datasetIndex) {
        int numPoints, degree;
        inputFile >> numPoints >> degree;

        points.clear();
        for (int i = 0; i < numPoints; ++i) {
            double x, y;
            inputFile >> x >> y;
            points.emplace_back(x, y);
        }

        BestSolution.second = 0;

        for (int generation = 0; generation < MaxGeneration; ++generation) {
            vector<pair<vector<double>, double>> sol_fitness;
            for (int i = 0; i < population; ++i) {
                vector<double> solution = GenerateRandomSolutions(degree);
                double fit = fitness(solution, degree);
                sol_fitness.emplace_back(solution, fit);
            }

            tournamentSelection(sol_fitness);
            CrossOver(offsprings, degree);
            Mutation(offsprings, degree, mutationLocation, generation);

            for (int i = 0; i < offsprings.size(); ++i) {
                double fit = fitness(offsprings[i], degree);
                sol_fitness.emplace_back(offsprings[i], fit);
            }

            tournamentSelection(sol_fitness);


            if (sol_fitness[0].second > BestSolution.second) {
                BestSolution = sol_fitness[0];
            }
        }


        outputFile << "Dataset: " << datasetIndex + 1 << endl;
        outputFile << "Coefficients: ";
        for (const auto &coeff : BestSolution.first) {
            outputFile << coeff << " ";
        }
        outputFile << endl;
        outputFile << "Mean Square Error: " << BestSolution.second << endl;
        outputFile << "---------------------------------" << endl;
    }

    inputFile.close();
    outputFile.close();

    return 0;
}
