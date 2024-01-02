#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "gmp.h"
#include "base64.h"
#include "buffer.h"
#include "bits.h"
#include "random.h"
#include "operating_modes.h"
#include "aes.h"
#include "version.h"
#include "network.h"
#include "dh.h"
#include "channel.h"
#include "certificate.h"

#define DEFAULT_SERVER_PORT 31415
#define DEFAULT_PORT 1789
#define DEFAULT_HOST "localhost"
#define DEFAULT_NAME "silviu-andrei.maftei"

static char *client_host;
static int client_port;
static char *user_name;

void handle_reply(char **from, int *portfrom, char **reply, char **packet){
    *packet = network_recv(1);
    if(!parse_packet(from, portfrom, reply, *packet)){
        return;
    };
    printf("Received \"%s\" from %s:%d!\n", *reply, *from, *portfrom);
}

void try_send(const char *host, const int port){
    char *hello = malloc(sizeof(char)*(40+strlen(user_name)+strlen(client_host)+1));
    sprintf(hello, "Hello! My name is %s, calling from %s.", user_name, client_host);
    network_send(host, port, client_host, client_port, hello);
    free(hello);
    char *packet = network_recv(1);
    if (packet == NULL) {
        printf("[ERROR] %s didn't reply to me (%s)\n", host, client_host);
        return;
    }
    free(packet); // discard first reply
    network_send(host, port, client_host, client_port, "I am the client sending to the server.\n");
}

void try_aes(){
    uchar *msg = (uchar*)"It's a long way to Tipperary";
    buffer_t clear, encrypted, key, IV, decrypted;
    mpz_t gab;

    buffer_init(&clear, strlen((char*)msg));
    buffer_init(&encrypted, 1);
    buffer_init(&key, BLOCK_LENGTH);
    buffer_init(&IV, BLOCK_LENGTH);
	
    mpz_init_set_str(gab, "12345612345678907890", 10);
    AES128_key_from_number(&key, gab);
    buffer_random(&IV, BLOCK_LENGTH);
    buffer_from_string(&clear, msg, strlen((char*)msg));

    aes_CBC_encrypt(&encrypted, &clear, &key, &IV, 's');

    buffer_init(&decrypted, 1);
    aes_CBC_decrypt(&decrypted, &encrypted, &key, 's');
    buffer_print(stdout, &decrypted);
    printf("\n");
	
    buffer_clear(&clear);
    buffer_clear(&encrypted);
    buffer_clear(&decrypted);
    buffer_clear(&key);
    buffer_clear(&IV);
    mpz_clear(gab);
}

void send_with_aes(const char *host, const int port, uchar *msg, mpz_t gab)
{
    //initialise utensils
    buffer_t plain_text, encrypted, InitValue, key, enc_in_64;
    const char *string_aux = "AES";
    buffer_init(&plain_text, strlen((char*)msg));
    buffer_init(&encrypted, 1);
    buffer_init(&InitValue, BLOCK_LENGTH);
    buffer_init(&key, BLOCK_LENGTH);
    buffer_init(&enc_in_64, 1);
    buffer_random(&InitValue, BLOCK_LENGTH);
    
    buffer_from_string(&plain_text, msg, strlen((char*)msg));
    
    //encrypt
    AES128_key_from_number(&key, gab);
    aes_CBC_encrypt(&encrypted, &plain_text, &key, &InitValue, 's'); 
    buffer_to_base64(&enc_in_64, &encrypted);
    
    //send to network
    network_send(host, port, client_host, client_port, string_aux);
    network_send(host, port, client_host, client_port, (char *)enc_in_64.tab);

    //tidy up
    buffer_clear(&plain_text);
    buffer_clear(&encrypted);
    buffer_clear(&InitValue);
    buffer_clear(&key);
    buffer_clear(&enc_in_64);
}

void try_send_aes(const char *host, const int port){
    uchar *msg = (uchar*)"It's a long way to Tipperary";
    mpz_t gab;

    mpz_init_set_str(gab, "12345612345678907890", 10);
    send_with_aes(host, port, msg, gab);
    mpz_clear(gab);
}

void prepare_cipher(buffer_t *encrypted, buffer_t *clear, buffer_t *key){
    buffer_t IV;
    buffer_init(&IV, BLOCK_LENGTH);
    buffer_random(&IV, BLOCK_LENGTH);
    aes_CBC_encrypt(encrypted, clear, key, &IV, 's');
    buffer_clear(&IV);
}

void encrypt_aes(buffer_t *encrypted, uchar *msg, mpz_t gab)
{
    //inits
    buffer_t plaintext, IV, key, aux;
    buffer_init(&plaintext, strlen((char*)msg));
    buffer_init(&IV, BLOCK_LENGTH);
    buffer_init(&key, BLOCK_LENGTH);
    buffer_init(&aux, 1);
    buffer_random(&IV, BLOCK_LENGTH);
    buffer_from_string(&plaintext, msg, strlen((char*)msg));

    //encrypt
    AES128_key_from_number(&key, gab);
    aes_CBC_encrypt(&aux, &plaintext, &key, &IV, 's'); 
    buffer_to_base64(encrypted, &aux);

    buffer_clear(&plaintext);
    buffer_clear(&IV);
    buffer_clear(&key);
    buffer_clear(&aux);
}

void CaseDH(const char *server_host, const int server_port, gmp_randstate_t state)
{
    //inits
    mpz_t p, g, a, ga, gb, gab;
    size_t nbits;
    char buf[1024], *msg, *packet;
    const char *string_aux = "AES";
    mpz_inits(p, g, a, ga, gb, gab, NULL);
    buffer_t aux;
    buffer_init(&aux, 1);
    channel_init(p, g);
    nbits = mpz_sizeinbase(p, 2)-1;
    DH_init(a, state, nbits);
    
    //step 1
    mpz_powm_sec(ga, g, a, p);
    msg_export_mpz(buf, "DH: ALICE/BOB CONNECT1 ", ga, 0);
    network_send(server_host, server_port, client_host, client_port, buf);
    packet = network_recv(-1);//listen indefinitely
    
    //step 2
    parse_packet(NULL, NULL, &msg, packet);
    if (msg_import_mpz(gb, msg, "DH: BOB/ALICE CONNECT2 ", 0) <= 0) return; // if message different break out of the function

    //step 3
    mpz_powm_sec(gab, gb, a, p);
    encrypt_aes(&aux, (uchar *)string_aux, gab);
    msg_export_string(buf, "DH: ALICE/BOB CONNECT3 ", (const char *)aux.tab);
    network_send(server_host, server_port, client_host, client_port, buf);
    
    //clean up
    mpz_clears(p, g, a, ga, gb, gab, NULL);
    buffer_clear(&aux);
    free(msg);
    free(packet);
}

int CaseSTS(const char *server_host, const int server_port, certificate_t *CA, mpz_t NA, mpz_t dA, mpz_t N_aut, mpz_t e_aut, gmp_randstate_t state)
{
    int retno=0;
    mpz_t p, g, a, ga, gb, gab, sA, sB, auxB;
    size_t nbits;
    char buf[1024], *tmp, *msg, *packet;
    buffer_t y, auxbuf, sigB, plaintext, enc, IV;
    certificate_t CB;
    int msg_received = 0;

    mpz_inits(p, g, ga, a, NULL);

    channel_init(p, g);
    nbits = mpz_sizeinbase(p, 2)-1;

    DH_init(a, state, nbits);
    mpz_powm_sec(ga, g, a, p);
    msg_export_mpz(buf, "STS: ALICE/BOB CONNECT1 ", ga, 0);
    network_send(server_host, server_port, client_host, client_port, buf);

    packet = network_recv(-1);
    parse_packet(NULL, NULL, &msg, packet);
    msg_received = msg_import_string(buf, msg, "STS: BOB/ALICE CONNECT2 ");
    if (msg_received <= 0) 
    {
        free(packet);
        free(msg);
        mpz_clears(p, g, ga, a, NULL);
        return 0;
    }
    free(packet);
    free(msg);

    buffer_init(&auxbuf, 1);
    buffer_init(&y, 1);
    buffer_from_string(&auxbuf, (uchar *)buf, strlen(buf));
    buffer_from_base64(&y, &auxbuf);
    buffer_clear(&auxbuf);

    packet = network_recv(-1);
    parse_packet(NULL, NULL, &msg, packet);
    mpz_inits(gb, gab, NULL);
    msg_received = msg_import_mpz(gb, msg, "STS: BOB/ALICE CONNECT2 ", 0);
    if (msg_received <= 0) 
    {
        free(packet);
        free(msg);
        buffer_clear(&y);
        mpz_clears(p, g, ga, a, gb, gab, NULL);
        return 0;
    }
    mpz_powm_sec(gab, gb, a, p);
    free(packet);
    free(msg);

    packet = network_recv(-1);
    parse_packet(NULL, NULL, &msg, packet);
    msg_received = msg_import_string(buf, msg, "STS: BOB/ALICE CONNECT2 ");
    if (msg_received <= 0) 
    {
        free(packet);
        free(msg);
        buffer_clear(&y);
        mpz_clears(p, g, ga, a, gb, gab, NULL);
        return 0;
    }
    free(packet);
    free(msg);
    init_certificate(&CB);
    certificate_from_string(&CB, buf);
    if (!valid_certificate(&CB, N_aut, e_aut)) 
    {
        printf("certificate invalid\n");
        mpz_clears(p, g, ga, a, gb, gab, NULL);
        buffer_clear(&y);
        clear_certificate(&CB);
        return retno;
    }
	
    buffer_init(&auxbuf, 1);
    buffer_init(&sigB, 1);
    mpz_inits(sB, auxB, NULL);
    buffer_from_mpz(&auxbuf, gab);
    aes_CBC_decrypt(&sigB, &y, &auxbuf, 's');
    buffer_to_mpz(sB, &sigB);
    concatenate_gb_ga(auxB, ga, gb, p);
    mpz_powm_sec(sB, sB, CB.e, CB.N);
    if (!mpz_cmp(auxB, sB)) 
    {
        mpz_clears(p, g, ga, a, gb, gab, sB, auxB, NULL);
        buffer_clear(&y);
        buffer_clear(&auxbuf);
        buffer_clear(&sigB);
        clear_certificate(&CB);
        return retno;
    }
    mpz_clears(sB, auxB, NULL);
    buffer_clear(&sigB);

    mpz_init(sA);
    buffer_init(&IV, BLOCK_LENGTH);
    buffer_init(&enc, 1);
    buffer_init(&plaintext, 1);
    SIGNSK(sA, ga, gb, p, NA, dA);
    buffer_random(&IV, BLOCK_LENGTH);
    buffer_from_mpz(&plaintext, sA);
    aes_CBC_encrypt(&enc, &plaintext, &auxbuf, &IV, 's');
    buffer_to_base64(&plaintext, &enc);
    tmp = (char *)string_from_buffer(&plaintext);
    mpz_clear(sA);
    buffer_clear(&IV);

    msg_export_string(buf, "STS: ALICE/BOB CONNECT3 ", tmp);
    free(tmp);

    network_send(server_host, server_port, client_host, client_port, buf);

    tmp = (char*)string_from_certificate(CA);
    msg_export_string(buf, "STS: ALICE/BOB CONNECT3 ", tmp);
    free(tmp);
    network_send(server_host, server_port, client_host, client_port, buf);

	
    mpz_clears(p, g, ga, a, gb, gab, NULL);
    buffer_clear(&y);
    clear_certificate(&CB);
	
    packet = network_recv(-1);
    parse_packet(NULL, NULL, &msg, packet);
    msg_received = msg_import_string(buf, msg, "STS: BOB/ALICE CONNECT4 ");
    if (msg_received <= 0) 
    {
        free(packet);
        free(msg);
        buffer_clear(&enc);
        buffer_clear(&plaintext);
        buffer_clear(&auxbuf);
        return retno;
    }
    buffer_from_string(&plaintext, (uchar *)buf, strlen(buf));
    buffer_from_base64(&enc, &plaintext);
    aes_CBC_decrypt(&plaintext, &enc, &auxbuf, 's');
    printf("last message: %s\n", plaintext.tab);
    
    free(packet);
    free(msg);
    buffer_clear(&enc);
    buffer_clear(&plaintext);
    buffer_clear(&auxbuf);
    return 1;
}

void CaptureTheFlag(const char *server_host, const int server_port,
                    certificate_t *CA, mpz_t NA, mpz_t dA, mpz_t N_aut, mpz_t e_aut,
                    gmp_randstate_t state){
    // Ask to capture the flag
    char ctf[] = "CTF: CONNECT";
    network_send(server_host, server_port, client_host, client_port, ctf);
    network_send(server_host, server_port, client_host, client_port, user_name);
    if (!CaseSTS(server_host, server_port, CA, NA, dA, N_aut, e_aut, state)) {
        printf("[CTF] ERROR, try again.\n");
        return;
    }
    char *packet = network_recv(5);
    char *from, *msg;
    parse_packet(&from, NULL, &msg, packet);
    if (strcmp(from, server_host) != 0) {
        printf("[CTF] You've been hacked by %s!\n", from);
        free(packet);
        free(from);
        free(msg);
        return;
    }
    mpz_t secret;
    mpz_init(secret);
    printf("Message = %s\n", msg);
    if (mpz_set_str(secret, msg, 16) == 0) {
        gmp_printf("[CTF] Congratulations ! You captured your flag!\nSecret=%#Zx\n", secret);
    }
    else {
        printf("[CTF] ERROR, try again.\n");
    }
    free(packet);
    free(from);
    free(msg);
    mpz_clear(secret);
}


void Usage(char *s){
    fprintf(stderr, "Client for version %s\n\n", VERSION);
    fprintf(stderr, "Usage: %s\t[--sendto SERVER_HOST (default %s)] [--port SERVER_PORT (default %d)]\n", s, DEFAULT_HOST, DEFAULT_SERVER_PORT);
    /* fprintf(stderr, ""); */
    fprintf(stderr, "\t\t[--hostname CLIENT_HOST (default %s)]\n", DEFAULT_HOST);
    fprintf(stderr, "\t\t[--listen CLIENT_PORT (default %d)]\n", DEFAULT_PORT);
    fprintf(stderr, "\t\t[ try_send | try_aes | try_send_aes | try_DH | try_STS | try_CTF ]\n");
    fprintf(stderr, "\t\t[--help]");
    fprintf(stderr, " [--name NAME]\n");
    fprintf(stderr, "\t\t[ OPTIONAL FILES ]\n");

    fprintf(stderr, "\nArguments:\n");
    fprintf(stderr, "\tNAME: \t\tOf the form name.lastname, lowercase [Example: %s].\n", DEFAULT_NAME);
    fprintf(stderr, "\ttry_send: \tTry to send a basic message. Modify and play with this function.\n");
    fprintf(stderr, "\ttry_aes: \tTry to encrypt a basic message with AES. Don't send anything.\n");
    fprintf(stderr, "\ttry_send_aes: \tEncrypt and send a basic message with AES.\n");
    fprintf(stderr, "\ttry_DH: \tPerform DH key exchange and encrypt a message with the shared key.\n");
    fprintf(stderr, "\ttry_STS: \tPerform Station-To-Station protocol with the server.\n");
    fprintf(stderr, "\ttry_CTF: \tTry to retrieve the secret flag prepared for you. Needs your name.\n");

    fprintf(stderr, "\nOptional files (Required for STS and CTF):\n\n");
    fprintf(stderr, "\tclient_certificate.txt\n\tclient_sk.txt\n\tauth_pk.txt\n");
}

int main(int argc, char *argv[])
{
    int opt=0;
    char *server_host;
    int server_host_defined = 0, server_port_defined=0;
    int client_host_defined=0, client_port_defined=0;
    int user_name_defined=0;
    int sts_flag=0, capture_the_flag=0;
    int server_port;
    FILE *in;
    certificate_t CA;
    mpz_t NA, dA, N_aut, e_aut;
    char *packet, *reply, *from;
    int portfrom;
    gmp_randstate_t state;

    static int random_port_flag;

    static struct option long_options[] = {
        {"random", no_argument, &random_port_flag, 1},
        {"help", no_argument, 0, 'h'},
        {"sendto", required_argument, 0, 's'},
        {"port", required_argument, 0, 'p'},
        {"hostname", required_argument, 0, 'c'},
        {"listen", required_argument, 0, 'l'},
        {"name", required_argument, 0, 'n'},
        {NULL, 0, 0, '\0'}
    };

    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "hs:p:c:l:n:",
                              long_options, &long_index )) != -1) {
        switch (opt) {
            case 0:
                break;
            case 'h':
                Usage(argv[0]);
                return 0;
            case 's':
                server_host = malloc(sizeof(char)*strlen(optarg)+1);
                strcpy(server_host, optarg);
                server_host_defined = 1;
                break;
            case 'p':
                server_port = atoi(optarg);
                server_port_defined = 1;
                break;
            case 'c':
                client_host = malloc(sizeof(char)*strlen(optarg)+1);
                strcpy(client_host, optarg);
                client_host_defined = 1;
                break;
            case 'l':
                client_port = atoi(optarg);
                client_port_defined = 1;
                break;
            case 'n':
                user_name = malloc(sizeof(char)*strlen(optarg)+1);
                strcpy(user_name, optarg);
                user_name_defined=1;
                break;
            case '?':
                fprintf(stderr, "Try %s --help\n", argv[0]);
                return 1;
            default:
                Usage(argv[0]);
                return 0;
        }

    }

    if (!server_host_defined) {
        server_host = malloc(sizeof(char)*(strlen(DEFAULT_HOST))+1);
        sprintf(server_host, "%s", DEFAULT_HOST);
        server_host_defined=1;
    }

    if (!server_port_defined) {
        server_port = DEFAULT_SERVER_PORT;
        server_port_defined=1;
    }

    if (!server_host_defined || !server_port_defined ) {
        Usage(argv[0]);
        return 0;
    }
    if (!client_host_defined) {
        client_host = malloc(sizeof(char)*(strlen(DEFAULT_HOST))+1);
        sprintf(client_host, "%s", DEFAULT_HOST);
    }
    if (!client_port_defined) {
        client_port = DEFAULT_PORT;
    }

    if (!user_name_defined) {
        user_name = malloc(sizeof(char)*(strlen(DEFAULT_NAME))+1);
        sprintf(user_name, "%s", DEFAULT_NAME);
        user_name_defined=0;
    }

    if (random_port_flag || check_port(client_port) == 0) {
        client_port_defined=1;
        srand(random_seed());
        do {
            client_port = (rand() % (65535 - 1024 + 1) + 1024);
        } while (check_port(client_port) == 0);
    }

    network_init(client_port);
    mpz_inits(NA, dA, N_aut, e_aut, NULL);
    gmp_randinit_default(state);

    if (argc == optind) {
        Usage(argv[0]);
        goto clear;
        return 0;
    }
    else {
        printf("Client for version %s\n", VERSION);
        if (strcmp(argv[optind], "try_send") == 0){
            try_send(server_host, server_port);
            handle_reply(&from, &portfrom, &reply, &packet);
            free(from);
            free(reply);
            free(packet);
            goto clear;
        }
        if(strcmp(argv[optind], "try_send_aes") == 0){
            try_send_aes(server_host, server_port);
            handle_reply(&from, &portfrom, &reply, &packet);
            free(from);
            free(reply);
            free(packet);
            goto clear;
        }
        if(strcmp(argv[optind], "try_aes") == 0){
            try_aes();
            goto clear;
        }
        if(strcmp(argv[optind], "try_DH") == 0){
            printf("Go for DH: \n");
            CaseDH(server_host, server_port, state);
            char *packet2 = network_recv(5);
            if (packet2 != NULL) {
                char *reply2;
                parse_packet(NULL, NULL, &reply2, packet2);
                if (strcmp(reply2, "DH: OK")==0) {
                    printf("try_DH: [OK]\n");

                }
                else {
                    printf("try_DH: [FAILED]\n");
                }
                free(packet2);
                free(reply2);
            }
            else {
                printf("try_DH: [FAILED]\n");
            }
            goto clear;
        }


        if(strcmp(argv[optind], "try_STS") == 0) {
            sts_flag=1;
        }
        if(strcmp(argv[optind], "try_CTF") == 0){
            capture_the_flag=1;
        }
        if (sts_flag || capture_the_flag) {
            if (argc-optind<3) {
                fprintf(stderr, "\nWrong arguments. Try %s --help\n\n", argv[0]);
                goto clear;
            }
            if((in = fopen(argv[++optind], "r")) == NULL){ // client_certificate.txt
                perror(argv[optind]);
                return -1;
            }
            optind++;
            init_certificate(&CA);
            extract_certificate(&CA, in);
            fclose(in);
#if DEBUG > 0
            print_certificate(&CA);
#endif

            if((in = fopen(argv[optind], "r")) == NULL){ // client_secret_key
                perror(argv[optind]);
                return -1;
            }
            optind++;
            read_secret_keys(NA, dA, in);
#if DEBUG > 0
            gmp_printf("\nNA=%Zd,\ndA=%Zd\n\n", NA, dA);
#endif
            fclose(in);

            if((in = fopen(argv[optind], "r")) == NULL){ // auth_public_key
                perror(argv[optind]);
                return -1;
            }
            optind++;
            read_public_keys(N_aut, e_aut, in);
            fclose(in);
#if DEBUG > 0
            gmp_printf("\nN_aut:=%Zd;\ne_aut:=%Zd;\n\n", N_aut, e_aut);
#endif
            if (!valid_certificate(&CA, N_aut, e_aut)) {
                fprintf(stderr, "[ERROR] Client certificate is not valid!\n");
                exit(-1);
            }
#if DEBUG > 0
            printf("valid certificate for CA-client? %d\n",
                   valid_certificate(&CA, N_aut, e_aut));
#endif
            if (sts_flag) {
                printf("Go for STS: \n");
                CaseSTS(server_host, server_port, &CA, NA, dA, N_aut, e_aut, state);
            }
            else if (capture_the_flag) {
                CaptureTheFlag(server_host, server_port, &CA, NA, dA, N_aut, e_aut, state);
            }
            clear_certificate(&CA);
        }
    }
    clear:
       free(server_host);
       free(client_host);
       free(user_name);
       mpz_clears(NA, dA, N_aut, e_aut, NULL);
       gmp_randclear(state);
       network_clear();
       return 0;
}