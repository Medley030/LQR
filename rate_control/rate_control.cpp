/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file RateControl.cpp
 */

#include "rate_control.hpp"
#include <px4_platform_common/defines.h>
#include <iostream>

using namespace matrix;
using namespace std;

void RateControl::setGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_gain_p = P;
	_gain_i = I;
	_gain_d = D;
}

void RateControl::setSaturationStatus(const Vector<bool, 3> &saturation_positive,
				      const Vector<bool, 3> &saturation_negative)
{
	_control_allocator_saturation_positive = saturation_positive;
	_control_allocator_saturation_negative = saturation_negative;
}

Vector3f RateControl::update(const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     const float dt, const bool landed)
{
	// angular rates error
	Vector3f rate_error = rate_sp - rate;

	// PID control with feed forward
	const Vector3f torque = _gain_p.emult(rate_error) + _rate_int - _gain_d.emult(angular_accel) + _gain_ff.emult(rate_sp);

	// update integral only if we are not landed
	if (!landed) {
		updateIntegral(rate_error, dt);
	}

	return torque;
}


Vector3f RateControl::mcupdate(const Vector3f &att, const Vector3f &att_sp, const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     const float dt, const bool landed)
{
	// state feedback gain
	float _dataKa[18] = {10, 0, 0, 5.6945, 0, 0, 0, 10, 0, 0, 5.6945, 0, 0, 0, 2.4495, 0, 0, 1.8082};
	Matrix<float, 3, 6> Ka(_dataKa);

	Vector3f _att_error = att - att_sp;
	//cout << "roll error : " << _att_error(0) << ", pitch error : " << _att_error(1) << ", yaw error : " << _att_error(2) << '\n';
	//cout << "setpoint : " << att_sp(0) << ", " << att_sp(1) << ", " << att_sp(2) << '\n';
	//cout << "attitude" << att(0) << ", " << att(1) << ", " << att(2) << '\n';
	//cout << "roll error : " << _att_error(0) << ", pitch error : " << _att_error(1) << ", yaw error : " << _att_error(2) << '\n';
	//Vector3f new_rate_sp = Vector3f(0.1f*rate_sp(0), 0.1f*rate_sp(1), 0.1f*rate_sp(2));
	Vector3f _rate_error = rate - rate_sp;
	//cout << "from rate : " << rate_sp(0) << ", " << rate_sp(1) << ", " << rate_sp(2) << '\n';
	//cout << "roll error : " << _rate_error(0) << ", pitch error : " << _rate_error(1) << ", yaw error : " << _rate_error(2) << '\n';


	float _datazero[6] = {0, 0, 0, 0, 0, 0};
	Matrix<float, 6, 1> ea(_datazero);

	for(int i = 0; i < 3; i++){
		ea(i,0) = _att_error(i);
	}

	for(int i = 0; i < 3; i++){
		ea(i + 3,0) = _rate_error(i);
	}

	Matrix<float, 3, 1> Ua = -Ka*ea;


	const Vector3f torque = Vector3f(Ua(0,0), Ua(1,0), Ua(2,0));

	// update integral only if we are not landed
	// if (!landed) {
	// 	updateIntegral(rate_error, dt);
	// }
	//cout << "torque" << torque(0) << ", " << torque(1) << ", " << torque(2) << '\n';
	return 0.01f*torque;
}

void RateControl::updateIntegral(Vector3f &rate_error, const float dt)
{
	for (int i = 0; i < 3; i++) {
		// prevent further positive control saturation
		if (_control_allocator_saturation_positive(i)) {
			rate_error(i) = math::min(rate_error(i), 0.f);
		}

		// prevent further negative control saturation
		if (_control_allocator_saturation_negative(i)) {
			rate_error(i) = math::max(rate_error(i), 0.f);
		}

		// I term factor: reduce the I gain with increasing rate error.
		// This counteracts a non-linear effect where the integral builds up quickly upon a large setpoint
		// change (noticeable in a bounce-back effect after a flip).
		// The formula leads to a gradual decrease w/o steps, while only affecting the cases where it should:
		// with the parameter set to 400 degrees, up to 100 deg rate error, i_factor is almost 1 (having no effect),
		// and up to 200 deg error leads to <25% reduction of I.
		float i_factor = rate_error(i) / math::radians(400.f);
		i_factor = math::max(0.0f, 1.f - i_factor * i_factor);

		// Perform the integration using a first order method
		float rate_i = _rate_int(i) + i_factor * _gain_i(i) * rate_error(i) * dt;

		// do not propagate the result if out of range or invalid
		if (PX4_ISFINITE(rate_i)) {
			_rate_int(i) = math::constrain(rate_i, -_lim_int(i), _lim_int(i));
		}
	}
}

void RateControl::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
	rate_ctrl_status.rollspeed_integ = _rate_int(0);
	rate_ctrl_status.pitchspeed_integ = _rate_int(1);
	rate_ctrl_status.yawspeed_integ = _rate_int(2);
}
