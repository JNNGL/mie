#pragma once

#include <functional>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <queue>

namespace mie {

    typedef std::function<void()> Task;

    class ThreadPool {
    public:
        explicit ThreadPool(size_t threads = std::thread::hardware_concurrency());
        ~ThreadPool();

        void enqueue(Task task);
        void join();

    private:
        std::vector<std::thread> workThreads;
        std::queue<Task> taskQueue;
        std::mutex queueMutex;
        std::condition_variable cv;
        bool stop = false;
    };

}