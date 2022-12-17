# Vincent Hodgins
# EECE 560
# Final Project

# I chose to solve this project in Python rather than TCL as I wanted to write all of
# the functions used myself to gain a better understanding and I'm no good at TCL. 

## Functions contained are:

## MR_find_m(m, s=0) ##
# Factors out all twos from prime m, and returns m with all of its prime factors of two removed
# If s is not null, it will return tuple (s, k), where s is the number of 2's and k is 
# m with its prime factors of two removed.


## miller_rabin(n, override_trialcount=0, trace=0, prelim=1) ##
# Tests whether a number n is prime by the miller-rabin algorithm.
# global TRIALCOUNT is used for the number of trials to preform if 
# override_trialcount is not specified. 
# Trace is by defauly 0
# Setting trace to 1 prints out generally what is going on in the function as it computes
# Setting trace to 2 prints out everything that is going on. 
# These trace conventions carry on for the rest of the functions in this library

## test_mr(n) ## 
# Tests the accuracy of the miller rabin implementation for primes less than n
# was just used during debug, but could be used to demonstrate MR's effectiveness

## find_n_bit_prime(n, trace=0) ## 
# Checks numbers randomly in the n-bit range for primality. Returns the first prime
# it confidently finds.

## is_safe_prime(n, trace=0) ##
# Checks whether given n is a safe prime. Returns true or false

## find_safe_prime(n, trace=0) ## 
# Continously calls find_n_bit_prime(n), and when a prime is found, 
# applies is_safe_prime(n) to it to check if it is a safe prime.
# Returns the first safe prime found

## extended_euclidean_algorithm(a,b) ## 
# Applies the extended euclidean algorithm to a and b.
# returns a dict containing:
# key ["Bez"] : tuple of the two found Bezout coefficients
# key ["GCD"] : the greatest common denominator of a and b
# key ["QTS"] : tuple the quotients of a and b through EEA 

## modular_inverse(a,b) ## 
# Determines a's inverse mod b using the EEA

## rsa_genkeys(n, en_key=0, trace=0) ##
# Generates a dict containing public private key pair tuples from
# given bit length n, and optional encryption key en_key
# Returns a dict containing:
# key ["publickey"] : tuple of (encrypytion_key, modulus)
# key ["privatekey"]: tuple of (decrypytion_key, modulus) 

## string_to_charInt(M) ##
# Converts a string M to an integer representing M with specifications:
# A -> 10, Z -> 35
# a -> 36, z -> 61
# any other char -> 99
# This is used to encrypt strings, by making them into numbers

## charInt_to_string(M) ## 
# Converts an integer as defined in the previous function to a string with 
# 10 -> A, 35 -> Z
# 36 -> a, 61 -> z
# 99 -> " "
# Does this by taking every two digits in the int and converting each to a letter

## rsa_encrypt(M, pubkey) ## 
# Encrypts a string or int M by the given pubkey tuple returned from rsa_genkeys

## rsa_decrypt(C, privkey, words=0) ##
# Decrypts ciphertext C by privkey tuple
# If words set to other than null/0, translate the decrypted result
# from int to english under charInt_to_string before returning

## full_loop_cryptosystem(M,n=50,en_key=0, words=0, trace=0) ## 
# generates an n bit public/private keypair, with an optional arg en_key
# prints that pair to console
# encrypts a string or int M under that keypair and prints result to console
# decrypts that message and prints to console
# If optional arg trace=1, prints generally what is happening to the console
# trace=2 prints every single step to console
# this problem effectively solves questions 1 and 2 of the assignment all in one go

#                Example output of this function
             
#                 >>> full_loop_cryptosystem("Hello World") 
#                 Public Key:  (3, 920084918827348126785033087253)
#                 Private Key:  (613389945884897472172663707059, 920084918827348126785033087253)
#                 Encrypted Message: 271200579262711308338412339733
#                 Now decrypting: Hello  World


## generator_test(n,p) ##
# Tests if n is a generator under prime p
# Was going to make an automated search function, but I quickly found a generator 
# without brute forcing in my assignment, so I did not make one

## Q_res(n,p) ## 
#  returns wether or not n is a quadratic residue of p

## tonelli_shanks(n,p) ##
# Calculates the modular sqare root of n under modulus p by the tonelli_shanks alg


 
## This function was taken from stackoverflow and is not mine, but it is not really used in anything other
## than to test the efficiency of my miller rabin function for low primes
def SieveOfEratosthenes(num):
    plist = []
    prime = [True for i in range(num+1)]
# boolean array
    p = 2
    while (p * p <= num):
  
        # If prime[p] is not
        # changed, then it is a prime
        if (prime[p] == True):
  
            # Updating all multiples of p
            for i in range(p * p, num+1, p):
                prime[i] = False
        p += 1
  
    # Print all prime numbers
    for p in range(2, num+1):
        if prime[p]:
            plist.append(p)
    return plist  

low_primes = SieveOfEratosthenes(100)


## Miller - Rabin Primality Test

import random, time
random.seed(time.time())
TRIALCOUNT = 20                                  # Number of random a value trials to try -- ended up not using this global

def MR_find_m(n, s=0):
    k=1
    current_guess = (n-1)/2                     # Start factoring out 2's from n
    while (True):    
        if current_guess%2!=0:                  # If no more 2's in n 
            if s:
                return (k, int(current_guess))
            else:
                return int(current_guess)           # return (n-1)/(2**k) and k
        else:
            current_guess/=2                    # factor a two out of n                             
            k+=1                                # increment how many 2's found


def miller_rabin(n, override_trialcount=0, trace=0, prelim=1):
    ## Input Validation
    if (n<4):
        return True#raise Exception("Chosen n is too low")
    if (n%2==0):                                
        return False#raise Exception("Yeah, did you think that number was really prime..")

    ## Import globals 
    global TRIALCOUNT, low_primes

    if prelim:                                 # this ended up not saving any real amount of time. but i assume it does         
        for i in low_primes:                   # for large large values? but maybe not due to how expmod works
            if n%i==0:
                return False

    if not override_trialcount:                 # if we gave a keyword arg we use that, else default to global (unused)
        override_trialcount = TRIALCOUNT

    ## Step 1 : Find m
    m = MR_find_m(n)                            # determine our m factor 
    if trace==2:
        print("Found m as",m)

    ## Step 2 : Complete Trials
    for i in range(0,override_trialcount):      # repeat following trialcount times
        a = random.randint(2,n-2)               # choose a random a in range {2,n-2}
        if trace==2:
            print(f"Testing with random a={a}")
        #x = (a**m) % n
        x = pow(a,m,n)                          
        if x in [-1,1]:                         # oddity 1
            return True
        else:
            x = pow(a,2*m,n)                      
            if x==-1:                           # oddity 2 
                return True


    return False



## Accuracy check Miller_Rabin
def test_mr(n):
    plist = SieveOfEratosthenes(n)
    hits = 0; misses=0
    for i in plist:
        if miller_rabin(i):
            hits+=1
        elif not miller_rabin(i):
            misses+=1

    print("Hit rate:",hits/len(plist), " Miss rate:", misses/len(plist))


## Below executes a test to determine accuracy out of curiousity
#test_mr(20000)


### PROBLEM 1: 

## Generate a public/private key pair for the RSA cryptosystem. Do this by creating:

# A Safe 50 bit prime p and a 50 bit prime q

def find_n_bit_prime(n, trace=0):
    while (True):
        i = random.randint(2**(n-1)+1, (2**n) -1)      # Start with an n bit number                                       
        if trace==2:                                   
            print("Currently at : ", i)             
        if miller_rabin(i, override_trialcount=1, trace=(trace==2)):
            if trace==2:
                print("Possibly found one... double checking")
            if miller_rabin(i, override_trialcount=5, prelim=0):
                if trace==2:
                    print("Higher probability... triple checking")
                if miller_rabin(i, override_trialcount=50, prelim=0):
                    if trace==2:
                        print("Found prime.")
                        print(i)
                return i


def is_safe_prime(n, trace=0):             # Preform miller_rabin on (n-1)/2 to prove safety
    if miller_rabin(int((n-1)/2), override_trialcount=25, trace=trace):     # function defined above that returns bool
        if trace:               
            print(f"{n} is a safe prime, composed of prime:{int((n-1)/2)}")
        return True
    return False


def find_safe_prime(n, trace=0):
    while(True):
        p = find_n_bit_prime(n, trace=trace)
        if trace==2:
            print("\n               ------ Safety Check -------")
        if is_safe_prime(p, trace=trace):
            return p


# Run find_safe_prime(50) to find a safe prime of length 50 bits. 
# Set keyword argument trace=1 for more detail, and 2 for even more detail. 


## 1. - c. -- To find decryption key a modular inverse must be found: 
# Algorithm implementation instructions from https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode

def extended_euclidean_algorithm(a,b):
    old_r, r = a, b                                 # This is mostly copied directly from the wikipedia pseudocode, 
    old_s, s = 1, 0                                 # but in doing so i learned that python supports parallel assignments
    old_t, t = 0, 1                                 # pretty cool! 
    while (r!=0):
        quotient = old_r // r
        old_r, r = r, old_r - (quotient*r)
        old_s, s = s, old_s - (quotient*s)
        old_t, t = t, old_t - (quotient*t)

    return {                                        # Index this dict to grab desired result 
        "Bez" : (old_s, old_t),
        "GCD" : old_r,
        "QTS" : (t,s)
    }


def modular_inverse(a,b):
    bez = extended_euclidean_algorithm(a,b)['Bez']
    return bez[0] % b


# Combine the above steps to generate a public private rsa key pair 
def rsa_genkeys(n, en_key=0, trace=0):
    random.seed(time.time())
    primes = (find_safe_prime(n, trace=trace), find_safe_prime(n, trace=trace))
    modulus = primes[0] * primes[1]
    totient = (primes[0]-1) * (primes[1]-1)
    if trace:
        print(f"p={primes[0]}, q={primes[1]}\nn={modulus}, phi={totient}")
    if en_key==0:
        en_key = low_primes[random.randint(0,len(low_primes))]
    
    de_key = modular_inverse(en_key, totient)

    return {
        'publickey' : (en_key, modulus),
        'privatekey' : (de_key, modulus)
    }


# If a string was given instead of a number, this function puts it into an appropriate form
# by converting each letter into two digits in a number, preserving case
def string_to_charInt(M):
    s = ''
    for i in M:
        if i.isupper():
            s+= str(ord(i)-55)
        elif i.islower():
            s+= str(ord(i)-61)
        else:
            s+= '99'
    return int(s)

# this function preforms the opposite of the above 
def charInt_to_string(M):
    split_strings = []
    M = str(M)
    for i in range(0, len(M), 2):
        two_chars = M[i:i+2]
        split_strings.append(int(two_chars))

    s =''
    for item in split_strings:
        if item==99:
            s+= ' '
        if item<36:
            s+= str( chr(item+55))
        else:
            s+= str( chr(item+61))
    return s
    
# Encrypts message given message M and tuple publickey pubkey
def rsa_encrypt(M, pubkey):
    if not M.isnumeric():
        if not M.isnumeric():
            M = string_to_charInt(M)
    return pow(M, pubkey[0], pubkey[1])

# Decrypts message given message M and tuple privkey, along with keyword 
# arg words, set to 1 if expecting a string result rather than numerical
def rsa_decrypt(C, privkey, words=0):
    pt = pow(C, privkey[0], privkey[1])
    if words:
        return charInt_to_string(pt)
    else:
        return pt


# Preforms a loop of all of the previous steps, sufficient to answer questions 1 and 2
def full_loop_cryptosystem(M,n=50,en_key=0, trace=0):
    words = not M.isnumeric()
    keys = rsa_genkeys(n, en_key, trace=trace)
    print(f"Public Key:  {keys['publickey']}\nPrivate Key:  {keys['privatekey']}")
    ciphertext = rsa_encrypt(M, keys['publickey'])
    plaintext = rsa_decrypt(ciphertext, keys["privatekey"], words=words)
    print(f"Encrypted Message: {ciphertext}\nNow decrypting: {plaintext}")




## 3 Test generator

def generator_test(n,p):
     p1 = p-1
     q = int(p1/2)
     a = pow(n,2,p)
     b = pow(n,q,p)
     #print(f"{n}**2={a}\n{n}**q={b}")
     if (a!=1) and (b!=1):
        return True
     return False

## 5 Tonelli-Shanks Algorithm

def Q_res(n,p):                         # Quadratic Residue finder
    res = pow(n, int((p-1)/2), p)
    if res==1:
        return 1
    elif res==p-1:
        return -1
    else:
        return 0 
    

# Following wikipedia algorithm https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#The_algorithm
def tonelli_shanks(n,p):
    if Q_res(n,p)!=1:
        print("No roots exist."); return False
    
    # Factor out powers of two first
    s,q = MR_find_m(p, s=True)


    # Search for a quadratic nonresidue
    z = 0
    for i in range(2, p):
        if Q_res(i,p)==-1:
            z = i; break

    
    # Initialize Variables before loop
    m = s
    c = pow(z, q,p)
    t = pow(n, q,p)
    r = pow(n, int((q+1)/2),p)

    # Loop ## this part I just could not figure out for some reason so its from stackoverflow
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r  

