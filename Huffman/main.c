//
//  main.c
//  Huffman
//
//  Created by Locke Phan on 4/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

#include <CoreFoundation/CoreFoundation.h>

void * allocate ( size_t size ) {
    void *space = malloc(size);
    if ( space == NULL ) {
        fprintf(stderr, "Error allocating memory.\n");
        exit(EXIT_FAILURE);
    }
    return space;
}

char *code_word[256];

struct Heap_Node {
    size_t              count;  
    int                 character_number;
    struct Heap_Node    *left;
    struct Heap_Node    *right;
    struct Heap_Node    *parent;
};

typedef struct Heap_Node Heap_Node;

struct Stack_Node {
    char                value;
    struct Stack_Node   *previous;
    struct Stack_Node   *next;
};

typedef struct Stack_Node Stack_Node;

struct custom_stack_t {
    struct Stack_Node   *head;
    struct Stack_Node   *current;
    size_t              count;
};



void stack_push ( struct custom_stack_t *stack , char c ) {
    Stack_Node *newNode = (Stack_Node *) allocate(sizeof(Stack_Node));
    
    newNode->previous   = NULL;
    newNode->next       = NULL;
    newNode->value      = c;
    
    if (stack->head == NULL) {
        stack->head     = newNode;
        stack->current  = newNode;
    } else {
        stack->current->next    = newNode;
        newNode->previous       = stack->current;
        stack->current          = newNode;
    }
    stack->count++;
}

void stack_pop ( struct custom_stack_t *stack ) {
    if (stack->current == NULL) return;    
    
    Stack_Node *popped_Node = stack->current;
    
    if (stack->current == stack->head) {
        stack->current  = NULL;
        stack->head     = NULL;
    } else {
        stack->current          = stack->current->previous;
        stack->current->next    = NULL;
    }
    free(popped_Node);
    popped_Node = NULL; // just to be safe
    stack->count--;
}

char * stack_print ( struct custom_stack_t *stack ) {
    char *string = (char *) allocate(stack->count + 1);
    
    Stack_Node *p = stack->head;
    for(int i = 0 ; p != NULL ; p = p->next, i++)
        string[i] = p->value;
    
    string[stack->count] = '\0';
    return string;
}

void walk ( struct custom_stack_t *stack, Heap_Node *root ) {
    if ( root == NULL ) return;
    if ( root->character_number >= 0 ) {
        char *myString = stack_print(stack);
        code_word[root->character_number] = myString;
    } else {
        stack_push(stack, '0');
        walk(stack, root->left);
        stack_pop(stack);
        
        stack_push(stack, '1');
        walk(stack, root->right);
        stack_pop(stack);
    }
}


CFComparisonResult compare ( const void *ptr1, const void *ptr2, void *info ) {
    Heap_Node *node1 = (Heap_Node *) ptr1;
    Heap_Node *node2 = (Heap_Node *) ptr2;
    
    if ( node1->count > node2->count )
        return kCFCompareGreaterThan;
    else if ( node1->count < node2->count )
        return kCFCompareLessThan;
    else 
        return kCFCompareEqualTo;
}


void free_huffman_code_tree ( Heap_Node *node ) {
    Heap_Node *left = node->left;
    Heap_Node *right = node->right;
    free(node);
    
    if ( left != NULL )
        free_huffman_code_tree(left);
    if ( right != NULL )
        free_huffman_code_tree(right);
}

Heap_Node * create_huffman_code_tree ( int *counts ) {
    // Initializes binary heap to be used as a priority queue
    CFBinaryHeapCallBacks myCallBacks;
    myCallBacks.compare = compare;
    myCallBacks.retain = NULL;
    myCallBacks.release = NULL;
    CFBinaryHeapRef heap = CFBinaryHeapCreate(kCFAllocatorDefault, 0, &myCallBacks, NULL);
    
    // Loop over the array, create nodes, and
    // add them to the priority queue
    size_t total = 0;
    for (int i = 0 ; i < 256 ; i++) {
        if (counts[i] > 0) {
            total = total + counts[i];
            Heap_Node *heap_node        = (Heap_Node *) allocate(sizeof(Heap_Node));
            heap_node->count            = counts[i];
            heap_node->character_number = i;
            heap_node->left             = NULL;
            heap_node->right            = NULL;
            heap_node->parent           = NULL;
            
            CFBinaryHeapAddValue(heap, heap_node);
        }
    }
    
    // Run Huffman Code algorithm to construct a prefix-free code tree
    CFIndex queue_size = CFBinaryHeapGetCount(heap) - 1;
    for ( int i = 0 ; i < queue_size; i++ ) {
        Heap_Node *node = (Heap_Node *) allocate(sizeof(Heap_Node));
        
        Heap_Node *x = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
        CFBinaryHeapRemoveMinimumValue(heap);
        Heap_Node *y = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
        CFBinaryHeapRemoveMinimumValue(heap);
        
        node->left = x;
        x->parent = node;
        node->right = y;
        y->parent = node;
        node->character_number = -1;
        node->count = x->count + y->count;
        
        CFBinaryHeapAddValue(heap, node);
    }
    
    // Extract the huffman code tree from the heap
    Heap_Node *node = (Heap_Node *) CFBinaryHeapGetMinimum(heap);
    CFRelease(heap);
    
    return node;
}

void encode (unsigned char *data, unsigned char *encodedData, size_t ord_data_size) {
    unsigned char byte = 0;
    size_t result_position = 0;
    size_t bit_position = 7;
    
    for ( int i = 0 ; i < ord_data_size ; i++ ) {
        unsigned char c = data[i];
        
        char *code = code_word[c];
        
        size_t code_length = strlen(code);
        for ( int j = 0 ; j < code_length ; j++ ) {
            if ( code[j] == '1' ) 
                byte |= ( 1 << bit_position );
            if ( bit_position == 0 ) {
                encodedData[result_position] = byte;
                result_position++;
                bit_position = 7;
                byte = 0;
                continue;
            }
            bit_position--;
        }
    }
    if ( byte != 0 )
        encodedData[result_position] = byte;
    
}

void decode (Heap_Node *root, unsigned char *decoded_data, unsigned char *encoded_data, size_t maxbits ) {    
    unsigned char byte;
    Heap_Node *current = root;
    size_t totalbits = 0;
    size_t decode_data_index = 0;
    size_t num_encoded_bytes = maxbits / 8;
    
    for ( int i = 0 ; i < num_encoded_bytes ; i++ ) {
        byte = encoded_data[i];
        
        for ( int j = 7 ; j >= 0 && totalbits < maxbits ; j--, totalbits++ ) {
            if ( byte & ( 1 << j ))
                current = current->right;
            else
                current = current->left;
            if ( current->character_number >= 0 ) {
                decoded_data[decode_data_index] = current->character_number;
                current = root;
                decode_data_index++;
            }
        }
    }
}

ssize_t get_file_size ( char *path ) {
    struct stat buffer;
    int error;
    if ( (error = stat(path, &buffer)) < 0) {
        perror("Function: get_file_size.");
        exit(EXIT_FAILURE);
    }
    
    return buffer.st_size;
}

ssize_t read_file_to_buffer ( char *path, void *buffer, size_t filesize ) {
    int fd;
    if ((fd = open(path, O_RDONLY)) < 0) {
        perror("Read file");
        exit(EXIT_FAILURE);
    }
    
    ssize_t bytes_read;
    if ((bytes_read = read(fd, buffer, filesize)) < 0) {
        perror("Can't read.");
        exit(EXIT_FAILURE);
    }
    close(fd);
    return bytes_read;
}

int main (int argc, const char * argv[])
{
        char *filepath = "/Users/tienloc47/test.txt";
    //    char *filepath = "/Users/tienloc47/sample_large_text_file.txt";
    //    char *filepath = "/Users/tienloc47/movie.avi";
    //    char *filepath = "/Users/tienloc47/random_file";
    //    char *filepath = "/Volumes/Time Machine/Movies/An.Edu.DVDSCR.XviD.avi";
    //    char *filepath = "/Volumes/Time Machine/Movies/City.Of.Angels.1998.DVDRip.Xvid.AC3.Magpie.avi";
    //    char *filepath = "/Volumes/Time Machine/Movies/Dear.John.2010.720p.BluRay.AC3.x264-UNiT3D.mkv";
    //char *filepath = "/Volumes/Time Machine/Movies/Fired.Up.2009.720p.BluRay.AC3.x264-UNiT3D.mkv";
    
    ssize_t filesize = get_file_size(filepath);
    
    unsigned char *ORG_DATA = (unsigned char *) allocate(filesize);
    
    ssize_t max = SSIZE_MAX;
    printf("%ld \n", max);
    
    printf("Read %ld bytes into buffer.\n", read_file_to_buffer(filepath, ORG_DATA, filesize));
    
    
    // Initializes counts array to store the number of each character
    int *counts = (int *) allocate(sizeof(int)*256); 
    memset(counts, 0, sizeof(int)*256);
    
    
    // Counts the number of each char in file
    unsigned char c;
    for ( int i = 0; i < filesize ; i++ ) {
        c = ORG_DATA[i];
        ++counts[c];
    }
    
    Heap_Node *node = create_huffman_code_tree(counts);
    
    // Initialize a stack to keep track of where we are in the tree
    struct custom_stack_t *stack = (struct custom_stack_t *) allocate(sizeof(struct custom_stack_t));
    stack->head = NULL;
    stack->current = NULL;
    stack->count = 0;
    
    // Traverse through the tree to collect our codewords
    walk(stack, node);
    
    
    size_t total_compressed_numbits = 0;
    for (int i = 0; i < 256 ; i++ ) {
        char *string = code_word[i];
        if (string == NULL) continue;
        size_t len = strlen(string);
        size_t size = len * counts[i];
        total_compressed_numbits += size;
    }
    
    printf("Size of file: %lu bytes.\n", node->count);
    printf("Size after compression: %lu bytes.\n", total_compressed_numbits/8);
    
    size_t compressed_byte_size;
    if (total_compressed_numbits % 8 == 0) {
        compressed_byte_size = total_compressed_numbits / 8;
    } else {
        compressed_byte_size = (total_compressed_numbits / 8) + 1;
    }
    
    unsigned char *compressed_data = (unsigned char *) allocate(compressed_byte_size);
    memset(compressed_data, 0, compressed_byte_size);
    printf("Start encoding process...\n");
    encode(ORG_DATA, compressed_data, filesize);
    
    unsigned char *decoded_data = (unsigned char *) allocate(filesize);
    printf("Start decoding process...\n");
    decode(node, decoded_data, compressed_data, total_compressed_numbits);
    
    for ( int i = 0; i < filesize ; i++ ) {
        if ( ORG_DATA[i] != decoded_data[i] )
            printf("index: %d ORG: %d , decoded: %d \n", i, ORG_DATA[i], decoded_data[i]);
    }
    
    
    // Release encode array
    for ( int i = 0 ; i < 256 ; i++ )
        free(code_word[i]);
    
    free(counts);    
    free(stack);
    free_huffman_code_tree(node);
    
    printf("Program finished executing.\n");
    return 0;
}

