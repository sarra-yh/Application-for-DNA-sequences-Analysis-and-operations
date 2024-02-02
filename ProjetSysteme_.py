import random
import numpy as np
def lire_adn_de_fichier(chemin_fichier):
    try:
        with open(chemin_fichier, 'r') as fichier:
            return fichier.read()
    except FileNotFoundError:
        print("Le fichier spécifié n'a pas été trouvé.")
        return None

def generer_adn(longueur):
    bases_adn = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases_adn) for _ in range(longueur))

def verifier_adn(adn):
    bases_adn = ['A', 'C', 'G', 'T']
    b=True
    for i in range(len(adn)) :
        if adn[i] not in bases_adn:
              b=False
    return  b

def arn(adn):
    s2 = ""
    for base in adn:
        if base == "T":
            s2 += "U"
        else:
            s2 += base
    return s2

def traduire_en_proteines(arn):
    code_genetique = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    proteines = ""
    for i in range(0, len(arn), 3):
        codon = arn[i:i+3]
        if codon in code_genetique:
            proteines += code_genetique[codon]

    return proteines   

def comp_inverse(adn):
    s2 = ""
    k = len(adn)
    for i in range(k):
         if(adn[i] == "T"):
        
              s2 = s2+"A"
         if(adn[i] == "A"):
        
              s2 = s2+"T"
         if(adn[i] == "C"):
        
              s2 = s2+"G"
         if(adn[i] == "G"):
        
              s2 = s2+"C"

    
    return s2[::-1]
def taux_gc(adn):
    
    d=dict(( ('A', 0),
           ('C', 0), 
           ('G', 0),
           ('T', 0) 
            ))

    for s in adn:
        d[s] += 1

    n= (100*(d['C'] + d['G']))/len(adn)
    return n

def freq_codon(adn):
    occ = {}
    freq={}
    total_codon=0
    for i in range(0,len(adn),3):
        occurence=0
        codon=adn[i:i+3]
        #pour chaque codon on va parcourir la sequence pour voir le nombre d'occurence de son apparition
        for j in range(0,len(adn),3):
            if adn[j:j+3] == codon :
                occurence=occurence+1
        occ[codon]=occurence
        total_codon=total_codon+1 #calculer le nombre total des codons
    for i in range(0,len(adn),3): #calculer la frequence pour chaque codon freq=(occ(codon)/nbrtotal)*100
          codon=adn[i:i+3]
          freq[codon]=round((occ[codon]/total_codon)*100,2)   #utiliser la fonction round pour minimiser le float résultant trés long
  
    return freq

def random_substitution_mutation(adn, num_mutations):
    
    bases_adn = ['A', 'C', 'G', 'T']
    
    mutated_sequence = list(adn)

    for _ in range(num_mutations):
        index = random.randint(0, len(adn) - 1)
        current_base = adn[index]
        new_base = random.choice([base for base in bases_adn if base != current_base])
        mutated_sequence[index] = new_base

    return ''.join(mutated_sequence)
###############trouver les motifs ADN############
def find_substring_locations(dna_string, substring):
   
    locations = []
    start = 0

    while start < len(dna_string):
        index = dna_string.find(substring, start)
        if index == -1:
            break
        locations.append(index + 1)  # Adding 1 to convert from 0-based index to 1-based index
        start = index + 1

    return locations
#####################"ADN CONSENSUS ET matrice profile#############
def generate_profile_matrix_single_sequence(dna_sequence):
    
    length = len(dna_sequence)
    profile_matrix = np.zeros((4, length), dtype=int)

    for i, nucleotide in enumerate(dna_sequence):
        if nucleotide == 'A':
            profile_matrix[0, i] += 1
        elif nucleotide == 'C':
            profile_matrix[1, i] += 1
        elif nucleotide == 'G':
            profile_matrix[2, i] += 1
        elif nucleotide == 'T':
            profile_matrix[3, i] += 1

    consensus_string = ''
    for i in range(length):
        max_index = np.argmax(profile_matrix[:, i])
        consensus_string += 'ACGT'[max_index]

    return consensus_string, profile_matrix
############sauvegarder les resultats

def sauvegarder_resultats(resultats, nom_fichier):
    try:
        with open(nom_fichier, 'w') as fichier:
            fichier.write(resultats)
        print(f"Résultats sauvegardés dans le fichier '{nom_fichier}'.")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde : {e}")


def main():
   
    
    
    choix_adn = int(input("voulez vous un utiliser un fichier ou generer aleatoirement votre sequence d'ADN?"
                      "1-fichier"
                      "2-ADN aléatoire"))
    if choix_adn == 1:
        chemin_fichier = input("Entrez le chemin du fichier contenant la chaîne ADN : ")
        adn_fichier = lire_adn_de_fichier(chemin_fichier)
        adn=adn_fichier
    elif choix_adn == 2 :
         longueur = int(input("Entrez la longueur de la chaîne ADN aléatoire : "))
         adn_genere = generer_adn(longueur)
         adn=adn_genere
    else:
        print("veuillez choisir le bon nombre") 
    
    # Début du programme


    while True:
        choix_utilisateur = int(input("Choisissez une option :\n"
                                    "1. Vérifier la validité de la séquence ADN\n"
                                    "2. Transcrire une séquence ADN en ARN\n"
                                    "3. Construire des protéines à partir d'ARN\n"
                                    "4. Calculer le complément inverse de l'ADN\n"
                                    "5. Calculer le taux de GC dans l'ADN\n"
                                    "6. Calculer la fréquence des codons dans l'ADN\n"
                                    "7. mutation adn\n"
                                    "8. motif adn\n"
                                    "9. génerer la matrice profile et consensus\n"
                                    "10. sauvegarder dans un fichier\n"
                                    "13. quitter le programme\n"
                                    "Entrez le numéro de votre choix : "))

        # Effectuer l'action en fonction du choix de l'utilisateur
        if choix_utilisateur == 1:
            verifier = verifier_adn(adn)
            print(verifier)
            resultats = f"Résultat de la vérification de validité : {verifier}"
        elif choix_utilisateur == 2:
            print(arn(adn))
            resultats = f"Séquence ARN : {arn(adn)}"
        elif choix_utilisateur == 3:
            arn_exemple = arn(adn)
            proteines_resultantes = traduire_en_proteines(arn_exemple)
            print("Séquence protéique résultante :", proteines_resultantes)
            resultats = f"Séquence protéique résultante : {proteines_resultantes}"
        elif choix_utilisateur == 4:
            inver = comp_inverse(adn)
            print(inver)
            resultats = f"Complément inverse de l'ADN : {inver}"

        elif choix_utilisateur == 5:
            taux = taux_gc(adn)
            resultats = f"Taux de GC dans l'ADN : {taux}%"
            print(f"{taux}%")
        elif choix_utilisateur == 6:
            print(freq_codon(adn))
            freq = freq_codon(adn)
            resultats = f"Fréquence des codons dans l'ADN : {freq}"
        elif choix_utilisateur == 7:
            num_mutations = int(input("Entrez le nombre de mutations à réaliser : "))
            mutated_adn = random_substitution_mutation(adn, num_mutations)
            print(f"ADN muté : {mutated_adn}")
            resultats = f"ADN muté : {mutated_adn}"
        elif choix_utilisateur == 8:
            substring = input("Entrer la sequence Motif que vous voulez chercher")
            locations = find_substring_locations(adn, substring)
            print("Locations of the substring:", locations)
            resultats = f"Locations du motif : {locations}"
        elif choix_utilisateur == 9:
            consensus, profile = generate_profile_matrix_single_sequence(adn)
            print("Consensus string:", consensus)
            print("Profile matrix:")
            print(profile)
            resultats = f"Consensus string : {consensus}\nProfile matrix :\n{profile}"
        elif choix_utilisateur == 10:
            nom_fichier = input("Entrez le nom du fichier pour sauvegarder les résultats : ")
            sauvegarder_resultats(resultats, nom_fichier)
        elif choix_utilisateur == 13:
            print("Quitter le programme.")
            break
            
        else:
            print("Choix invalide. Veuillez entrer un numéro de choix valide.")

        
    
if __name__ == "__main__":
    main()
