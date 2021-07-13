#include "../common/allocator.h"
#include <pthread.h>
pthread_mutex_t locks = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t empty = PTHREAD_COND_INITIALIZER;
template <class T>
struct Node
{
  T value;
  Node<T>* next;
};

extern std::atomic<bool> no_more_enqueues;

template <class T>
class OneLockBlockingQueue
{
    Node<T>* q_head;
    Node<T>* q_tail;
    std::atomic<bool> wakeup_dq;
    CustomAllocator my_allocator_;
public:
    OneLockBlockingQueue() : my_allocator_()
    {
        std::cout << "Using OneLockBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next = NULL;
        //my_allocator_.freeNode(newNode);
        q_head = newNode;
        q_tail = newNode;
        wakeup_dq.store(false);
    }

    void enqueue(T value)
    {
      Node<T>* node = (Node<T>* )my_allocator_.newNode();
      node->value = value;
      node->next = NULL;
      //Append to q_tail and update the queue
      pthread_mutex_lock(&locks);
      q_tail->next = node;
      q_tail = node;
      wakeup_dq.store(true);
      pthread_mutex_unlock(&locks);
    }

    bool dequeue(T *value)
    {

        pthread_mutex_lock(&locks);
        Node<T>* node = q_head;
        Node<T>* new_head = node->next;

        while(new_head==NULL){
          // if(no_more_enqueues.load()==false){
          //   std::cout<<"Start waiting"<<std::endl;
          //   pthread_cond_wait(&empty,&locks);
          //   std::cout<<"Stop waiting"<<std::endl;
          // }
          wakeup_dq.store(false);
          pthread_mutex_unlock(&locks);
          while(no_more_enqueues.load()==false && wakeup_dq.load()==false){

          }
          pthread_mutex_lock(&locks);
          new_head = q_head->next;
          if(new_head!=NULL){
            node = q_head;
            wakeup_dq.store(true);
          }else if(no_more_enqueues.load()==true){
            pthread_mutex_unlock(&locks);
            return false;
          }

        }
        *value = new_head->value;
        q_head = new_head;
        pthread_mutex_unlock(&locks);
        my_allocator_.freeNode(node);
        return true;
    }

    void cleanup()
    {
      pthread_mutex_destroy(&locks);
      pthread_cond_destroy(&empty);
        my_allocator_.cleanup();
    }
};
