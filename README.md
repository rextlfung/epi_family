# Test 2D EPI fully sampled, and SMS slice profile

GE phantom

P,epi2d,20230912_UN_EPI,3.7:
git commit a7c0b0feffd8842b2511e80b8d76e570f132c438
with reference shots

## 16 Sep 2023

### epi 2d multislice scan, 10 slices, 3.5mm

P,epi2d,Nz=10,16Sep2023.7

- slice locations ~match scanner epi
- can't tell in this phantom if +freq produces offset in S or I

```
recon_epi2d.m
```


### 2D GRE image of SMS slice profile

Create gre2dmb.tar:
```
write2DGRE;
```

P,gre2dmb.7

To recon and display slice profile:
```
recon2DGRE;
```


