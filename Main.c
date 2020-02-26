#include <stdo.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// constantes d ela funcion y el algoritmo genetico 
#define POPULATION_SIZE 100; // creamos poblaciones
#define FX_LOWER_BOUND -20;
#define FX_UPPER_BOUND 20;
#define PRECISION 3;

// definicion de un individuo 
typedef struct ind_t {
    int* chromosome;
    double x;
    double fitness; // aptitud
    int parents[2]; // quienes son los padres de ese individuo 
    int mutation_place; // punto de mutacion para la cruza de los padres
    int crossover_place; 
} Individual; // esto es un individuo 

//poblaciones, ruleta y elititista 
Individual* parents;
Individual* offspring;
Individual* the_best; // elitista
double* roulette;

//variables utilies 
unsigned chromosome_length;
double crossover_probablity;
double mutation_probability;
unsigned max_generations;
unsigned selected_father; // -----> CORRECCION!!
unsigned selected_mother; // -----> CORRECCION!!

/*  +-------------------------------------------+
    |   PARAMETROS DEL ALGORITMO   |
    +-------------------------------------------+
*/

void getParametros(){
    printf("\nNumero maximo de generaciones: "); scanf("%u",&max_generations);
    printf("\nProbrabilidad de cruza: "); scanf("%lf",&crossover_probablity);
    printf("\nProbrabilidad de mutacion: "); scanf("%lf",&mutation_probability);
}

/*  +-------------------------------------------+
    |   RESERVAR MEMORIA PARA LAS POBLACIONES   |
    +-------------------------------------------+
*/
void allocateMemory() {
    unsigned required_byte = sizeof(Individual) * POPULATION_SIZE;
    parents = (Individual*) malloc(required_bytes);
    offspring = (Individual*) malloc(required_bytes);

    int i;
    chromosome_length = ceil(log2((FX_UPPER_BOUND - FX_LOWER_BOUND) * pow(10, PRECISION)));
    
    for(i=0; i < POPULATION_SIZE; i++) {
        // ---------> CORRECCION DE AQUI EN ADELANTE!!
        parents[i].chromosome = (int*) calloc(chromosome_length, sizeof(int));
        offspring[i].chromosome = (int*) calloc(chromosome_length, sizeof(int));
        parents[i].x = offspring[i].x = RAND_MAX;
        parents[i].fitness = offspring[i].fitness = 0;
        parents[i].parents[0] = offspring[i].parents[0] = - 1;
        parents[i].parents[1] = offspring[i].parents[1] = - 1;
        parents[i].crossover_place = offspring[i].crossover_place = -1;
        parents[i].mutation_place = offspring[i].mutation_place = -1;  
    }
    the_best.chomosome = (int*) calloc(chromosome_length, size(int));
    the_best.fitness = 0; // AGREGADO
    roulette = (double*) malloc(sizeof(double) * POPULATION_SIZE);
}
/*  +-------------------------------------------+
    |   DEVOLVER UN ALEATORIO SOBRE UN INTERVALO|
    +-------------------------------------------+
*/
double randomDouble(double a, double b) {
    return (b - a) * (ran() / (double)RAND_MAX) + a;
}

/*  +-------------------------------------------+
    |   SIMULAR LANZADO DE UNA MONEDA AL AIRE   |
    +-------------------------------------------+
*/
int flip(double probability) {
    int valor = 0;
    if(randomDouble(0,1) <= probability) { // ---------> CORRECCION!!!
        valor = 1;
    }
    return valor;
}

/*  +-------------------------------------------+
    |   INICIALIZAMOS LA PRIMERA GENERACION     |
    +-------------------------------------------+
*/
void createFristGeneration() {
     int i, j;
     for( i=0; i< POPULATION_SIZE; i++) {
         for(j=0; j < chromosome_length; j++){
             parents[i].chromosome[j] = flip(0.5);
         }
     }
}

/*  +--------------------------------------------------------------------------+
    |   DECODIFICAR EL GENOTIPO EN BINARIO A REAL Y APLICAR AJUSTE AL RANGO    |
    +--------------------------------------------------------------------------+
*/
double binary2real(int* chromosome) {
    int i;
    double aux = 0.0;
    for(i = chromosome_length - 1; i >= 0; i--) {
        if(chromosome[i] == 1) { // CORRECTION!!
            aux += (double) pow(2, chromosome_length - i - 1);
        }
    }
    return FX_LOWER_BOUND + ( (aux * (FX_UPPER_BOUND - FX_LOWER_BOUND)) / (pow(2, chromosome_length) - 1));
} 

/*  +--------------------------------------------------------------------------+
    |   DECODIFICAR EL GENOTIPO EN BINARIO A REAL Y APLICAR AJUSTE AL RANGO    |
    +--------------------------------------------------------------------------+
*/
void evaluateTargetFunction(Individual* individual) {
    individual->x = binary2real(individual->chromosome);
    individual->fitness = 1 / (pow(individual->x, 2) + 0.001);
}
//evaluar a la poblacion 
/*  +-------------------------------------------+
    |   EVALUAR APTITUDE DE UNA POBLACION       |
    +-------------------------------------------+
*/
void evaluatePopulation(Individual* population) {
    int i;
    for(i=0; i POPULATION_SIZE;i++) {
        evaluateTargetFunction(&population[i]);
    }
}
/*  +----------------------------------------------------------------------------------+
    |   LLENAR LA RULETA CON LA PROBABILIDAD DE CADA INVIVIDUO PARA SER SELECCIONADO   |
    +----------------------------------------------------------------------------------+
*/
void updateRoulette(){
    int i;
    double sum_fitness = 0.0;
    for(i = 0; i < POPULATION_SIZE; i++){
        sum_fitness += popilation[i].fitness;
    }
    for(i = 0; i < POPULATION_SIZE; i++){
        roulette[i] = population[i].fitness / sum_fitness;
    }
}

/*  +-------------------------------------+
    |   HACEMOS LA SELECCION POR RULETA   |
    +-------------------------------------+
*/
unsigned rouletteWheelSelection(){
    double r = randomDouble(0,1);
    double sum = 0.0;
    int i, current_individual; //---------> CORRECTION DESDE AQUI!!!!
    for(i = POPULATION_SIZE; sum < r; i++) {
        current_individual = i % POPULATION_SIZE;
        sum += roulette[current_individual];
    }
    return current_individual;
}
/*  +------------------------------------------------+
    |   RECOMBINACION DE LOS PADRES SELECCIONADOS    |
    +------------------------------------------------+
*/
void crossover(Individual* father, Individual* mother, Individual* child1, Individual* child2){
    int i;
    if(flip(crossover_probability)){
        unsigned p = (unsigned) randomDouble(1, chromosome_length - 2);
        for(i = 0; i <= p; i++) {
            child1->chromosome[i] = father->chromosome[i];
            child2->chromosome[ p + 1 ] = mother->chromosome[i];
        }
        for(i = p+1; i < chromosome_length; i++) {
            child1->chromosome[i] = mother->chromosome[i];
            child2->chromosome[ i - p - 1 ] = father->chromosome[i];
        }
        child1->crossover_place = child2->crossover_place = p;
        child1->parents[0] = child2->parents[0] = selected_father + 1; // CORRECCION !!
        child1->parents[1] = child2->parents[1] = selected_mother + 1; // COREECCION !!
    }
    else{
        for(i = 0; i < chromosome_length; i++) {
            child1->chromosome[i] = father->chromosome[i];
            child2->chromosome[ i ] = mother->chromosome[i];
        }
        child1->crossover_place = child2->crossover_place = - 1;
        child1->parents[0] = child2->parents[0] = 0; // CORRECCION !!
        child1->parents[1] = child2->parents[1] = 0; // COREECCION !!
    }
}
/*  +-------------------------------------------+
    |   REALIZAR MUTACION                       |
    +-------------------------------------------+
*/
void mutation(Individual* individual) { 
    if(flip(mutation_probability)) {
        unsigned p = (unsigned) randomDouble(0, chromosome_length - 1);
        individual->chromosome[p] = 1 - individual->chromosome[p];
        individual->mutacion_place = p;
    }
    else{
        individual->mutacion_place = -1;
    }
}

/*  +--------------------------------------------------------------------------+
    |   Ley del mas fuerte                                                     |
    +--------------------------------------------------------------------------+
*/
// problacion elitistas

void elitism() {
    unsigned best_parent; // CORRECCION !!
    unsigned worst_child1, worsd_child2;
    int i;
    best_parent = worst_child1 = worsd2 =0;

    for(i=0; i < POPULATION_SIZE; i++) {
        if( offspring[i].fitness < offspring[worsd_child1].fitness) {
            worsd_child1 =i;
        }
        else{
            if( offspring[i].fitness < offspring[worsd_child2].fitness) {
            worsd_child2 = i;
        }
        if(parents[i].fitness > parents[best_parent].fitness){
            best_parent = i;
        }
    }

    offspring[worst_child1] = parents[best_parent];
    offspring[worst_child2] = parents[best_parent];

}
/*  +--------------------------------------------------------------------------+
    |   IMPRIMIR CROMOSOMA CON INDICADORES                                     |
    +--------------------------------------------------------------------------+
*/
void printChromosome(Individual* individual){
    int i;
    for(i=0; i < chromosome_length; i++){
        if(i == individual->mutation_place) printf("(");
        printf("%d", indivifual->chromosome[i]);
        if(i == individual->mutation_place) printf(")");
        if(i == individual->crossover_place) printf("/");
    }
}
/*  +--------------------------------------------------------------------------+
    |   MOSTRAR LA INFORMACION DE UNA POBLACION                                |
    +--------------------------------------------------------------------------+
*/
void printPopulationDetail(Individual* population){

    int i, current_best=0; //-----> AGREGADO
    double fitness_avg = 0.0 //----> AGREGADO
    printf("\n\n----------------------------------------------------------------\n");
    printf("#\tChromosome\tx\tFitness\tParents")
    printf("\n----------------------------------------------------------------\n");
    for(i = 0; i < POPULATION_SIZE; i++){
        printf("\n%03d ", i+1);
        printChromosome(&population[i]);
        printf("%.3f\t%.3f\t(%d,%d)", population[i].x, population[i].fitness,
                                      population[i].parents[0], population[i].parents[1]);
    //---------------------------------AGREGADO ----------------------------------------
        fitness_avg += population[i].fitness;
        if(population[i].fitness > population[current_best].fitness) {
            current_best = i;
        }
        if(population[current_best].fitness > the_best.fitness) {
            the_best = population[i];
        }
    }
    fitness_arv /= POPULATION_SIZE;
    printf("\n----------------------------------------------------------------\n");
    printf("\nAptitud promedio (pobalcion): %.3f", fitness_avg);
    printf("\n Mejor aptitude: %.3f", population[current_best].fitness);
}

// Funcion principal 
int main() {
    //preparacion
    getParameters();
    srand((long)tiem(NULL));
    allocateMemory();
    createFristGeneration();
    evaluatePopulation(parents);
    Individual* temp_helper; // -------> CORRECION!!!!
    int generation, i;
    for(generation =0; generation < max_generation; generation++){
        updateRoulette(parents);
        printPopulationDetail(parents);
        // proceso generacion 
        for(i =0; i < POPULATION_SIZE-1; i+=2){
            selected_father = rouletteWheelSelection();
            selected_mother = rouletteWheelSelection();
            crossover(&parents[selected_father], &parents[selected_mother], & offspring[i], &offspring[i+1]);
            mutation(&offspring[i]);
            mutation(&offspring[i+1]);
            evaluateTargetFunction(&offspring[i]);
            evaluateTargetFunction(&offspring[i+1]);
        }
        elitism();
        temp_helper = parents;
        parents = offspring;
        offspring = tem_helper;
        printf("\n\n\tTermino con exito generation %d\n", generation+1);
    }
    //--------------------------------AGREGADO-----------------------------------------------
    printf("\n\n**************************************************************************");
    printf("\n\t\t\tEL MEJOR DE TODOS");
    printf("\n**************************************************************************");
    printf("\n\tCadena binaria: "); printChromosome(&the_best);
    printf("\n\t %.3f\tFitness = %.3f", the_best.x, the_best.fitness);
    printf("\n\tPadres:  (%d, %d)\n", the_best.parents[0], the_best.parents[1]);
    free(parents);
    free(offspring);
    free(roulette);
    return 0;
}