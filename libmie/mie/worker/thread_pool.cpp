#include "thread_pool.h"

namespace mie {

    ThreadPool::ThreadPool(size_t threads) {
        for (size_t i = 0; i < threads; i++) {
            workThreads.emplace_back([this] {
                while (true) {
                    Task task;

                    {
                        std::unique_lock lock(queueMutex);

                        cv.wait(lock, [this] {
                            return !taskQueue.empty() || stop;
                        });

                        if (taskQueue.empty()) {
                            if (stop) {
                                return;
                            }
                            continue;
                        }

                        task = std::move(taskQueue.front());
                        taskQueue.pop();
                    }

                    task();
                }
            });
        }
    }

    ThreadPool::~ThreadPool() {
        join();
    }

    void ThreadPool::join() {
        if (stop) {
            return;
        }

        {
            std::unique_lock lock(queueMutex);
            stop = true;
        }

        cv.notify_all();

        for (std::thread& thread : workThreads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }

    void ThreadPool::enqueue(Task task) {
        {
            std::unique_lock lock(queueMutex);
            taskQueue.emplace(std::move(task));
        }

        cv.notify_one();
    }

}