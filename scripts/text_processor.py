#!/usr/bin/env python3

import queue
import threading

import time

__author__ = 'Bence Cserna'

sentences = queue.Queue()
output = queue.Queue()


def worker():
    while True:
        item = sentences.get()
        if item is None:
            break

        output.put(item=process(item))
        sentences.task_done()


def process(sentence):
    return [s.split() for s in sentence]


def main():
    sentences_list = [[" cu ius affert adversarium. Eam suscipit phaedrum at. Qui ex semper verear inciderint." for _
                      in range(0, 100000)] for _ in range(0, 16)]

    # Start the worker threads
    threads = []
    worker_threads = 4

    for i in range(worker_threads):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    # Add all sentences to the queue
    for sentence in sentences_list:
        sentences.put(sentence)

    # block until all tasks are done
    sentences.join()

    # stop workers
    for i in range(worker_threads):
        sentences.put(None)
    for t in threads:
        t.join()

    print("Done")


if __name__ == "__main__":
    start = time.time()
    main()
    print("Execution time {}".format(time.time() - start))
