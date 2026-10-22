"""Small alternative event networks, adapted to the GAGG simulation benchmark.

The U-Net follows the multiscale encoder/decoder idea, without the published
LSTM bottleneck. The Transformer uses local token attention and dense heads,
not the original SDD paper's event encoding. Neither is a paper reproduction.
Both consume the same five measured-voltage features as the frozen CNN.
"""
from __future__ import annotations
import math
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F


class DoubleConv(nn.Module):
    def __init__(self, incoming, outgoing):
        super().__init__()
        self.layers=nn.Sequential(nn.Conv1d(incoming,outgoing,5,padding=2),nn.SiLU(),
                                  nn.Conv1d(outgoing,outgoing,5,padding=2),nn.SiLU())
    def forward(self,x):
        return self.layers(x)


class UNetEvent(nn.Module):
    """Four downsample levels and skip connections; three dense event heads."""
    def __init__(self):
        super().__init__()
        sizes=(12,20,32,48,64)
        self.encoders=nn.ModuleList([DoubleConv(5 if i==0 else sizes[i-1],c)
                                     for i,c in enumerate(sizes)])
        self.decoders=nn.ModuleList([DoubleConv(sizes[i]+sizes[i-1],sizes[i-1])
                                     for i in range(len(sizes)-1,0,-1)])
        self.head=nn.Conv1d(sizes[0],3,1)
    def forward(self,x):
        skips=[]
        for i,encoder in enumerate(self.encoders):
            if i:x=F.avg_pool1d(x,2)
            x=encoder(x);skips.append(x)
        for decoder,skip in zip(self.decoders,reversed(skips[:-1])):
            x=F.interpolate(x,size=skip.shape[-1],mode='linear',align_corners=False)
            x=decoder(torch.cat((x,skip),dim=1))
        return self.head(x)


class TransformerEvent(nn.Module):
    """8 ns tokens, 32 embedding dimensions, four heads and two encoder layers.

    A native-grid convolution branch retains fine local features. Inference
    uses the same 1024-sample local windows as training, avoiding quadratic
    attention on 10000 samples and extrapolation of positional coordinates.
    """
    def __init__(self):
        super().__init__()
        self.embedding=nn.Conv1d(5,32,9,stride=4,padding=4)
        layer=nn.TransformerEncoderLayer(d_model=32,nhead=4,dim_feedforward=96,
            dropout=0.,activation='gelu',batch_first=True,norm_first=True)
        self.encoder=nn.TransformerEncoder(layer,2,enable_nested_tensor=False)
        # TransformerEncoder clones the prototype; initialize independent weights.
        for parameter in self.encoder.parameters():
            if parameter.ndim>1:nn.init.xavier_uniform_(parameter)
        self.local=nn.Sequential(nn.Conv1d(5,8,5,padding=2),nn.SiLU())
        self.decoder=nn.Sequential(nn.Conv1d(40,24,5,padding=2),nn.SiLU(),
                                   nn.Conv1d(24,24,5,padding=2),nn.SiLU(),nn.Conv1d(24,3,1))
        position=torch.arange(256,dtype=torch.float32)[:,None]
        scale=torch.exp(torch.arange(0,32,2,dtype=torch.float32)*(-math.log(10000.)/32.))
        encoding=torch.empty(256,32)
        encoding[:,0::2]=torch.sin(position*scale);encoding[:,1::2]=torch.cos(position*scale)
        self.register_buffer('position',encoding[None],persistent=True)
    def forward(self,x):
        if x.shape[-1]!=1024:
            raise ValueError('Transformer inference must use its 1024-sample training window')
        tokens=self.embedding(x).transpose(1,2)
        encoded=self.encoder(tokens+self.position).transpose(1,2)
        encoded=F.interpolate(encoded,size=x.shape[-1],mode='linear',align_corners=False)
        return self.decoder(torch.cat((encoded,self.local(x)),dim=1))


def make_model(kind):
    if kind=='unet':return UNetEvent()
    if kind=='transformer':return TransformerEvent()
    raise ValueError(f'Unknown new architecture: {kind}')


def windowed_prediction(model,features,batch_size=8):
    """Predict only the guarded 500-point core of every 1024-point window.

    Each input point is emitted exactly once; no seam averaging or labels/rate
    inputs. The complete-record feature transform has already retained history.
    Zero padding is used only beyond the acquired record; scoring excludes edges.
    """
    if features.ndim!=2 or features.shape[0]!=5 or not np.isfinite(features).all():
        raise ValueError('Expected five finite measured-input feature channels')
    n=features.shape[-1];chunks=[];starts=list(range(0,n,500))
    padded=np.pad(features,((0,0),(262,762)),mode='constant')
    for start in starts:chunks.append(padded[:,start:start+1024])
    output=np.empty((3,n),dtype=np.float32)
    with torch.inference_mode():
        for first in range(0,len(chunks),batch_size):
            predicted=model(torch.from_numpy(np.stack(chunks[first:first+batch_size]))).cpu().numpy()
            for j,start in enumerate(starts[first:first+batch_size]):
                length=min(500,n-start)
                output[:,start:start+length]=predicted[j,:,262:262+length]
    return dict(probability=1/(1+np.exp(-np.clip(output[0],-30.,30.))),
                charge_pe=3000.*np.exp(np.clip(output[1],-5.,5.)),
                time_offset_ns=4.*np.tanh(output[2]))
