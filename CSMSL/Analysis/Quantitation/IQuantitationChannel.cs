﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CSMSL.Chemistry;

namespace CSMSL.Analysis.Quantitation
{
    public interface IQuantitationChannel: IMass
    {
        bool IsMS1Based { get; }     

    }
}