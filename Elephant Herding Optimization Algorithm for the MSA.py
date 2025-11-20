#EHO2 Algorithm for Multiple Sequence Alignment (MSA)
#Bio-inspired by Elephant Herding Optimization


import random
import numpy as np
import blosum
import time
import os



GAP_PENALTY = -2


class Elephant:
    def __init__(self, alignment):
        self.alignment = alignment  # List of aligned sequences
        self.fitness = None

    def __str__(self):
        return ''.join(self.alignment)

def read_fasta(filename):
    sequences = []
    with open(filename, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                    seq = ''
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences

def pad_sequences(sequences, pad_char='-'):
    max_len = max(len(seq) for seq in sequences)
    return [seq + pad_char * (max_len - len(seq)) for seq in sequences]

def calculate_fitness(alignment):
    global NFE
    try:
        NFE = NFE + 1
    except Exception:
        NFE = 1
    m = globals().get('matrixBlsm')
    if m is None:
        m = {"A": {"A": 4, "R": -1, "N": -2, "D": -2, "C": 0, "Q": -1, "E": -1, "G": 0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 0, "W": -3, "Y": -2, "V": 0}, "R": {"A": -1, "R": 5, "N": 0, "D": -2, "C": -3, "Q": 1, "E": 0, "G": -2, "H": 0, "I": -3, "L": -2, "K": 2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3}, "N": {"A": -2, "R": 0, "N": 6, "D": 1, "C": -3, "Q": 0, "E": 0, "G": 0, "H": 1, "I": -3, "L": -3, "K": 0, "M": -2, "F": -3, "P": -2, "S": 1, "T": 0, "W": -4, "Y": -2, "V": -3}, "D": {"A": -2, "R": -2, "N": 1, "D": 6, "C": -3, "Q": 0, "E": 2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S": 0, "T": -1, "W": -4, "Y": -3, "V": -3}, "C": {"A": 0, "R": -3, "N": -3, "D": -3, "C": 9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1}, "Q": {"A": -1, "R": 1, "N": 0, "D": 0, "C": -3, "Q": 5, "E": 2, "G": -2, "H": 0, "I": -3, "L": -2, "K": 1, "M": 0, "F": -3, "P": -1, "S": 0, "T": -1, "W": -2, "Y": -1, "V": -2}, "E": {"A": -1, "R": 0, "N": 0, "D": 2, "C": -4, "Q": 2, "E": 5, "G": -2, "H": 0, "I": -3, "L": -3, "K": 1, "M": -2, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2}, "G": {"A": 0, "R": -2, "N": 0, "D": -1, "C": -3, "Q": -2, "E": -2, "G": 6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0, "T": -2, "W": -2, "Y": -3, "V": -3}, "H": {"A": -2, "R": 0, "N": 1, "D": -1, "C": -3, "Q": 0, "E": 0, "G": -2, "H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2, "V": -3}, "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3}, "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I": 2, "L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1}, "K": {"A": -1, "R": 2, "N": 0, "D": -1, "C": -3, "Q": 1, "E": 1, "G": -2, "H": -1, "I": -3, "L": -2, "K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2}, "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q": 0, "E": -2, "G": -3, "H": -2, "I": 1, "L": 2, "K": -1, "M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1}, "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I": 0, "L": 0, "K": -3, "M": 0, "F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1}, "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2}, "S": {"A": 1, "R": -1, "N": 1, "D": 0, "C": -1, "Q": 0, "E": 0, "G": 0, "H": -1, "I": -2, "L": -2, "K": 0, "M": -1, "F": -2, "P": -1, "S": 4, "T": 1, "W": -3, "Y": -2, "V": -2}, "T": {"A": 0, "R": -1, "N": 0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 5, "W": -2, "Y": -2, "V": 0}, "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F": 1, "P": -4, "S": -3, "T": -2, "W": 11, "Y": 2, "V": -3}, "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H": 2, "I": -1, "L": -1, "K": -2, "M": -1, "F": 3, "P": -3, "S": -2, "T": -2, "W": 2, "Y": 7, "V": -1}, "V": {"A": 0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I": 3, "L": 1, "K": -2, "M": 1, "F": -1, "P": -2, "S": -2, "T": 0, "W": -3, "Y": -1, "V": 4}}
        globals()['matrixBlsm'] = m
    if not hasattr(calculate_fitness, 'min_score'):
        vals = []
        for a in m:
            for b in m[a]:
                v = m[a][b]
                if isinstance(v, (int, float)):
                    vals.append(v)
        calculate_fitness.min_score = min(vals) if vals else 0.0
        calculate_fitness.max_score = max(vals) if vals else 1.0
        if calculate_fitness.max_score == calculate_fitness.min_score:
            calculate_fitness.max_score = calculate_fitness.min_score + 1.0
    min_s = calculate_fitness.min_score
    max_s = calculate_fitness.max_score
    denom = (max_s - min_s)
    score = 0.0
    n = len(alignment)
    if n == 0:
        return 0.0
    L = len(alignment[0]) if alignment[0] else 0
    for i in range(n-1):
        row_i = alignment[i]
        for j in range(i+1, n):
            row_j = alignment[j]
            for k in range(L):
                a = row_i[k]
                b = row_j[k]
                if a != '-' and b != '-':
                    raw_score = globals()['matrixBlsm'][a][b]
                    norm = (raw_score - min_s) / denom
                    score = score + norm
    return score

def initialize_population(sequences, pop_size):
    population = []
    for _ in range(pop_size):
        alignment = pad_sequences(sequences)
        # Optionally shuffle columns for diversity
        alignment = [list(seq) for seq in alignment]
        for col in range(len(alignment[0])):
            if random.random() < 0.1:
                for row in alignment:
                    row[col] = '-'
        alignment = [''.join(row) for row in alignment]
        elephant = Elephant(alignment)
        elephant.fitness = calculate_fitness(alignment)
        population.append(elephant)
    return population

def group_into_clans(population, num_clans):
    population.sort(key=lambda e: e.fitness, reverse=True)
    clans = [[] for _ in range(num_clans)]
    for idx, elephant in enumerate(population):
        clans[idx % num_clans].append(elephant)
    return clans





def get_aligned_cut_index(seq, aa_cut_index):
    
    #Dada una secuencia alineada (con gaps) y un índice de corte de aminoácido,
    #devuelve el índice en la secuencia alineada donde ocurre ese corte.
    
    aa_count = 0
    for i, c in enumerate(seq):
        if c != '-':
            aa_count += 1
        if aa_count == aa_cut_index:
            return i + 1  # Corte después de este aminoácido
    return len(seq)  # Si el corte es al final



def clan_updating_operator(clan, crossProbabilty):
    
    #Operador de actualización del clan con crossover biológico.
    #Al final, iguala la longitud de todas las secuencias rellenando con gaps al final.
    
    if len(clan) <= 1:
        return
    matriarch = max(clan, key=lambda e: e.fitness)
    for elephant in clan:
        if elephant is not matriarch:
            # Realizar crossover biológico con probabilidad
            if random.random() < crossProbabilty:
                new_alignment = []
                for seq1, seq2 in zip(elephant.alignment, matriarch.alignment):
                    # Contar aminoácidos reales (sin gaps) en ambas secuencias
                    aa_seq1 = [c for c in seq1 if c != '-']
                    aa_seq2 = [c for c in seq2 if c != '-']
                    max_cut = min(len(aa_seq1), len(aa_seq2))
                    if max_cut <= 1:
                        new_alignment.append(seq1)  # No se puede cruzar
                        continue
                    # Elegir punto de corte biológico
                    cut_aa_index = random.randint(1, max_cut - 1)
                    # Mapear a índices alineados
                    cut1 = get_aligned_cut_index(seq1, cut_aa_index)
                    cut2 = get_aligned_cut_index(seq2, cut_aa_index)
                    # Hacer crossover en el mismo punto biológico
                    new_seq = seq1[:cut1] + seq2[cut2:]
                    new_alignment.append(new_seq)
                # Igualar longitudes rellenando con gaps al final
                max_len = max(len(seq) for seq in new_alignment)
                new_alignment = [seq + '-' * (max_len - len(seq)) for seq in new_alignment]
                elephant.alignment = new_alignment
                elephant.fitness = calculate_fitness(new_alignment)




def separation_operator(clan, sequences):   #Phase 2
   
    if len(clan) <= 1:
        return

    matriarch = max(clan, key=lambda e: e.fitness)
    bestFitness = matriarch.fitness



    # Crear nuevo alineamiento aleatorio a partir de las secuencias originales
    alignment = [list(seq) for seq in sequences]  # Convertir a listas para modificar
    newElephantFitness = 0


    for _ in range(20):  #ciclo de intentos de mejora

        # Insertar gaps aleatorios
        num_gaps_to_insert = random.randint(1, 3)  # Insertar entre 1 y 3 gaps

        for _ in range(num_gaps_to_insert):
            # Seleccionar una secuencia aleatoria
            seq_index = random.randint(0, len(alignment) - 1)
            # Seleccionar una posición aleatoria en esa secuencia
            pos = random.randint(0, len(alignment[seq_index]))
            # Insertar gap en esa posición (desplaza aminoácidos a la derecha)
            alignment[seq_index].insert(pos, '-')

        # Igualar longitudes: encontrar la secuencia más larga
        max_length = max(len(seq) for seq in alignment)

        # Agregar gaps al final de las secuencias más cortas
        for seq in alignment:
            while len(seq) < max_length:
                seq.append('-')

        # Convertir de vuelta a strings
        alignment = [''.join(seq) for seq in alignment]

        # Crear nuevo elefante con el alineamiento generado
        new_elephant = Elephant(alignment)
        new_elephant.fitness = calculate_fitness(alignment)
        newElephantFitness = new_elephant.fitness
        if newElephantFitness > bestFitness:        # Si el nuevo elefante mejora la fitness, lo aceptamos
            break
        else: #convierte las secuencias de new elephant a listas
            alignment = [list(seq) for seq in alignment]







    # Eliminar el peor elefante del clan
    worst = min(clan, key=lambda e: e.fitness)
    clan.remove(worst)
    clan.append(new_elephant)




#fasta_file = 'C:\secuenciasBFOA\multifasta.FASTA'  # Replace with your FASTA file
#fasta_file = 'IDE_protein_sequences.fasta'

def eho2_msa(fasta_file = 'yourFastaFile.fasta' ,scheme = "B", pop_size=36, num_clans=6, max_iter=88, crossProbabilty= 0.8):
    global time
    global original_sequences
    original_sequences = read_fasta(fasta_file)  # Read original sequences from FASTA

    sequences = read_fasta(fasta_file)
    population = initialize_population(sequences, pop_size)
    progress = []
    for iteration in range(max_iter):
        clans = group_into_clans(population, num_clans)  #en cada iteracion se reordenan los clanes
        for clan in clans:
            clan_updating_operator(clan, crossProbabilty)
        for clan in clans:
            separation_operator(clan, sequences)
        # Merge clans
        population = [elephant for clan in clans for elephant in clan]
        best_fitness = max(e.fitness for e in population)
        progress.append(best_fitness)
        #registra el tiempo
        elapsed_time = time.time() - start_time
        print('Iteration ,', iteration+1,',Sequences Set ,', fasta_file,', Scheme: ,', scheme, ', Iterations ,', max_iter, ', Probabilidad de Cruza ,' , crossProbabilty, ', Population ,', pop_size, ',Num Clans ,', num_clans,', Best fitness: ,', best_fitness, ', Cumulative NFE: , ' ,NFE ,', Cumulative time: ,' ,elapsed_time )
    best_elephant = max(population, key=lambda e: e.fitness)
    # print('Best alignment found:')
    # print(best_elephant)
    return best_elephant.alignment, progress




def remove_gaps(seq):
    
    # Elimina todos los gaps (-) de una secuencia alineada.
    
    return seq.replace('-', '')

def validate(original_sequences, final_alignment):

    # print("=== VALIDACIoN DE INTEGRIDAD BIOLoGICA ===")
    # print(f"Comparando {len(original_sequences)} secuencias originales vs alineadas...")
    # print()

    all_valid = True
    validation_results = []

    for i, (original, aligned) in enumerate(zip(original_sequences, final_alignment)):
        # Remover gaps de la secuencia alineada
        aligned_no_gaps = remove_gaps(aligned)

        # Comparar con la secuencia original
        is_valid = (original == aligned_no_gaps)

        validation_results.append({
            'sequence_id': i + 1,
            'original': original,
            'aligned': aligned,
            'aligned_no_gaps': aligned_no_gaps,
            'is_valid': is_valid,
            'original_length': len(original),
            'aligned_length': len(aligned),
            'gaps_count': aligned.count('-')
        })

        # if is_valid:
            # print(f"✓ Secuencia {i+1}: VaLIDA")
            # print(f"  Original:     {original}")
            # print(f"  Alineada:     {aligned}")
            # print(f"  Sin gaps:     {aligned_no_gaps}")
            # print(f"  Gaps anadidos: {aligned.count('-')}")
        # else:
        #     print(f"✗ Secuencia {i+1}: INVaLIDA - ¡INFORMACIoN PERDIDA!")
        #     print(f"  Original:     {original}")
        #     print(f"  Alineada:     {aligned}")
        #     print(f"  Sin gaps:     {aligned_no_gaps}")
        #     print(f"  Diferencia:   {set(original) - set(aligned_no_gaps)}")
        #     all_valid = False
        # print()

    # Resumen final
    # print("=== RESUMEN DE VALIDACIoN ===")
    if all_valid:
        print()
        # print("✓ VALIDACIoN EXITOSA: Toda la informacion biologica se conservo.")
        # print("✓ El alineamiento es biológicamente integro.")
    else:
        print("✗ ---------------------------            VALIDACIÓN FALLIDA: Se detecto perdida de informacion biologica.")
        # print("✗ El alineamiento NO es biologicamente integro.")

    # print(f"Total de secuencias validadas: {len(original_sequences)}")
    # print(f"Secuencias validas: {sum(1 for r in validation_results if r['is_valid'])}")
    # print(f"Secuencias invalidas: {sum(1 for r in validation_results if not r['is_valid'])}")

    return all_valid, validation_results


    return validate(original_sequences, best_elephant.alignment)[0]



if __name__ == '__main__':
    matrixBlsm = blosum.BLOSUM(62)
    NFE = 0  # Número de evaluaciones de la función objetivo
    #inicia la cuenta del tiempo
    start_time = time.time()






    best_elephant, progress = eho2_msa()

    # Ejemplo de validación al final del algoritmo (ajusta según tu flujo):
    # original_sequences = read_fasta('C:\secuenciasBFOA\multifasta.FASTA')
    # best_elephant = max(population, key=lambda e: e.fitness)  # Ajusta si tu variable es diferente
    is_valid, results = validate(original_sequences, best_elephant)
    if is_valid:
        print("El alineamiento final es biologicamente valido.")
    else:
        print("------------------------------        ¡Atencion! El alineamiento final tiene problemas de integridad.")