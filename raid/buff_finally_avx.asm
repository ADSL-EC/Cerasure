;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Copyright(c) 2011-2015 Intel Corporation All rights reserved.
;
;  Redistribution and use in source and binary forms, with or without
;  modification, are permitted provided that the following conditions
;  are met:
;    * Redistributions of source code must retain the above copyright
;      notice, this list of conditions and the following disclaimer.
;    * Redistributions in binary form must reproduce the above copyright
;      notice, this list of conditions and the following disclaimer in
;      the documentation and/or other materials provided with the
;      distribution.
;    * Neither the name of Intel Corporation nor the names of its
;      contributors may be used to endorse or promote products derived
;      from this software without specific prior written permission.
;
;  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
;  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
;  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
;  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
;  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Optimized xor of N source vectors using AVX
;;; int xor_gen_avx(int vects, int len, void **array)

;;; Generates xor parity vector from N (vects-1) sources in array of pointers
;;; (**array).  Last pointer is the dest.
;;; Vectors must be aligned to 32 bytes.  Length can be any value.

%include "reg_sizes.asm"

%ifidn __OUTPUT_FORMAT__, elf64
 %define arg0  rdi
 %define arg1  rsi
 %define arg2  rdx
;  %define arg3  rcx
;  %define arg4  r8
;  %define arg5  r9

 %define tmp   r11
 %define tmp2  r10
 %define tmp3  r12     ; must be saved and restored
 %define return rax
 %define func(x) x: endbranch
  %macro FUNC_SAVE 0
	push	r12
 %endmacro
 %macro FUNC_RESTORE 0
	pop	r12
 %endmacro


%elifidn __OUTPUT_FORMAT__, win64
 %define arg0  rcx
 %define arg1  rdx
 %define arg2  r8
 %define arg3  r9
 %define tmp   r11
 %define tmp2  r10
 %define tmp3  r13
 %define func(x) proc_frame x
 %define return rax
 %define stack_size  2*32 + 8 	;must be an odd multiple of 8

 %macro FUNC_SAVE 0
	alloc_stack	stack_size
	vmovdqu	[rsp + 0*32], ymm6
	vmovdqu	[rsp + 1*32], ymm7
	end_prolog
 %endmacro
 %macro FUNC_RESTORE 0
	vmovdqu	ymm6, [rsp + 0*32]
	vmovdqu	ymm7, [rsp + 1*32]
	add	rsp, stack_size
 %endmacro

%endif	;output formats


%define len arg0
%define offset arg1
%define src arg2

; %define len arg3
; %define src arg4
; %define dest arg5
%define pos tmp3
%define ptr tmp
%define vec_i tmp2
%define tmp2.b al
%define pos tmp3
%define PS 8

; %define vec arg0
; %define len arg1
; %define ptr arg3
; %define tmp2 rax
; %define tmp2.b al
; %define pos tmp3
; %define PS 8

;;; Use Non-temporal load/stor 
%ifdef NO_NT_LDST
 %define XLDR vmovdqa
 %define XSTR vmovdqa
%else
 %define XLDR vmovdqa
 %define XSTR vmovdqa
 %define XNSTR vmovntdq
%endif


default rel
[bits 64]

section .text

align 16
mk_global  buff_finally_avx, function
func(buff_finally_avx)
	FUNC_SAVE
    xor pos,pos
    XLDR ymm4,[offset]
    XLDR ymm5,[offset+32]
    XLDR ymm6,[offset+64]
    XLDR ymm7,[offset+96]
loop128:
    XLDR ymm0,[src+pos]
    XLDR ymm1,[src+pos+32]
    XLDR ymm2,[src+pos+64]
    XLDR ymm3,[src+pos+96]
      

    
    VPADDQ ymm0, ymm0, ymm4 
    VPADDQ ymm1, ymm1, ymm5 
    VPADDQ ymm2, ymm2, ymm6 
    VPADDQ ymm3, ymm3, ymm7 
    ; VPADDQ ymm4, ymm4, ymm9 
    ; VPADDQ ymm5, ymm5, ymm9 
    ; VPADDQ ymm6, ymm6, ymm9 
    ; VPADDQ ymm7, ymm7, ymm9 

    XSTR [src+pos],ymm0
    XSTR [src+pos+32],ymm1
    XSTR [src+pos+64],ymm2
    XSTR [src+pos+96],ymm3
    ; XNSTR [src+pos+128],ymm4
    ; XNSTR [src+pos+160],ymm5
    ; XNSTR [src+pos+192],ymm6
    ; XNSTR [src+pos+224],ymm7
    add	pos, 128
    cmp	pos, len
    jle	loop128
return_pass:
	FUNC_RESTORE
	mov	return, 0
	ret


return_fail:
	FUNC_RESTORE
	mov	return, 1
	ret

endproc_frame

section .data

;;;       func         core, ver, snum
slversion buff_finally_avx, 02,   05,  0037

