//
//  ExecutionService.h
//  ParallelBalls
//
//  Created by msmith on 9/9/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#ifndef __ParallelBalls__ExecutionService__
#define __ParallelBalls__ExecutionService__

#include <iostream>

#include <future>
#include <functional>
#include <thread>
#include <queue>
#include <vector>
#include <condition_variable>
#include <mutex>

class ExecutionThread{
private:
    std::condition_variable cv;
    std::mutex m;
    std::queue<std::function<void()>> tasks;
    bool running=false;
    std::thread* worker;
    void execution();

public:
    void start();
    void submit(std::function<void()> fun );
    
    ~ExecutionThread(){
    
    }
            
    
    
};

class ExecutionService{
    private:
        std::vector<ExecutionThread*> executors;
        std::queue<std::promise<int>*> promises;
        int next;
        int MAX;
    public:
    ExecutionService(int threads);
    void submit(const std::function<void()> &callable);
    void waitForExecution();
};

#endif /* defined(__ParallelBalls__ExecutionService__) */
