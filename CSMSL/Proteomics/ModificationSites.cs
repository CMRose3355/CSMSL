﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSMSL.Proteomics
{
    [Flags]
    public enum ModificationSites
    {
        None = 0,
        A = 1,
        R = 2,
        N = 4,
        D = 8,
        C = 16,
        E = 32,
        Q = 64,
        G = 128,
        H = 256,
        I = 512,
        L = 1024,
        K = 2048,
        M = 4096,
        F = 8192,
        P = 16384,
        S = 32768,
        T = 65536,
        U = 131072,
        W = 262144,
        Y = 524288,
        V = 1048576,
        NPep = 2097152,
        PepC = 4194304,
        NProt = 8388608,
        ProtC = 16777216,
        All = Int32.MaxValue,
    }
}