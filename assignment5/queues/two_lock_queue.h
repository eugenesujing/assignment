#include "../common/allocator.h"
#include <pthread.h>
pthread_mutex_t lock1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock2 = PTHREAD_MUTEX_INITIALIZER;

template <class T>
struct Node
{
  T value;
  Node<T>* next;
};

template <class T>
class TwoLockQueue
{
public:
  Node<T>* q_head;
  Node<T>* q_tail;
    CustomAllocator my_allocator_;

public:
    TwoLockQueue() : my_allocator_()
    {
        std::cout << "Using TwoLockQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        my_allocator_.freeNode(newNode);
        newNode->next = NULL;
        q_head = newNode;
        q_tail = newNode;

    }

    void enqueue(T value)
    {
      Node<T>* node = (Node<T>* )my_allocator_.newNode();
      node->value = value;
      node->next = NULL;
      //Append to q_tail and update the queue
      pthread_mutex_lock(&lock1);
      q_tail->next = node;
      q_tail = node;
      pthread_mutex_unlock(&lock1);
    }

    bool dequeue(T *value)
    {

        pthread_mutex_lock(&lock2);
        Node<T>* node = q_head;
        Node<T>* new_head = node->next;

        if(new_head==NULL){
          pthread_mutex_unlock(&lock2);
          return false;

        }
        *value = new_head->value;
        q_head = new_head;
        pthread_mutex_unlock(&lock2);
        my_allocator_.freeNode(node);
        return true;
    }
    void cleanup()
    {
      pthread_mutex_destroy(&lock1);
      pthread_mutex_destroy(&lock2);
        my_allocator_.cleanup();
    }
};
