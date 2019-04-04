#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:45:49 2019

@author: antoineleblevec
"""
# =============================================================================
# Classe pour compter le nombre de fichiers dans un dossier
# =============================================================================
class Counter(object):
    def __init__(self, path):
        if not path:
            raise ValueError('Its a trap')
        self.path = path
        self.files = 0

    def work(self):
        for entry in scandir(self.path):
            if entry.is_dir() and not entry.is_symlink():
                path = os.path.join(self.path, entry.name)
                counter = Counter(path)
                yield from counter.work()
            else:
                self.files += 1
        yield self

    def __str__(self):
       return '{} {}'.format(self.path, self.files)
   


# =============================================================================
# Compteur du nombres de fichiers
# =============================================================================
counter = f1.Counter(name_dir)
total = 0
for cls in counter.work():
    total += cls.files