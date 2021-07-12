#include "../common/allocator.h"
#include <pthread.h>

static pthread_mutex_t* locks = new pthread_mutex_t;

template <class T>
struct Node
{
  T value;
  Node<T>* next;
};

template <class T>
class OneLockQueue
{
    Node<T>* q_head;
    Node<T>* q_tail;
    CustomAllocator my_allocator_;
public:
    OneLockQueue() : my_allocator_()
    {
        std::cout << "Using OneLockQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
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
      pthread_mutex_lock(locks);
      q_tail->next = node;
      q_tail = node;
      pthread_mutex_unlock(locks);
    }

    bool dequeue(T *value)
    {
        bool ret_value = false;
        pthread_mutex_lock(locks);
        if(q_head->next!=NULL){
          Node<T>* node = q_head;
          *value = q_head->next->value;
          q_head = node->next;
          my_allocator_.freeNode(node);
          ret_value = true;
        }
        pthread_mutex_unlock(locks);
        return ret_value;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }
};
